cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.develop(path="$(homedir())/github/SimpleSDMLayers.jl")
using Revise
using SimpleSDMLayers
using Test

# Import QGIS-cropped layer as example
ref = geotiff(SimpleSDMPredictor, "../../rasters/wc1_cropped.tif")

## Step 2 : Fix getindex bounding coordinates

# Load things
layer = SimpleSDMPredictor(WorldClim, BioClim, 1)
layer_subset = layer[coords]
size(layer_subset) # should be 330, 570
left, right, bottom, top = coords

# Preliminary tests
layer_subset
ref
# Latitudes should be (20.083333333333332, 74.91666666666667)
# Longitudes should be (-144.91666666666666, -50.083333333333336)
isequal(latitudes(ref), latitudes(layer_subset)) # false
isapprox(latitudes(ref), latitudes(layer_subset)) # true

coords_subset = (layer_subset.left, layer_subset.right, layer_subset.bottom, layer_subset.top)
coords_ref = (ref.left, ref.right, ref.bottom, ref.top)
isequal.(coords_ref, coords_subset) # false
isapprox.(coords_ref, coords_subset) # true

## Investigate the bug
# 1. getindex(layer, coords)
imax = _match_longitude(layer, isnothing(right) ? layer.right : right; side=:right) # ok
imin = _match_longitude(layer, isnothing(left) ? layer.left : left; side=:left) # ok
jmax = _match_latitude(layer, isnothing(top) ? layer.top : top; side=:top) # ok
jmin = _match_latitude(layer, isnothing(bottom) ? layer.bottom : bottom; side=:bottom) # ok
layer[jmin:jmax, imin:imax] # problem is here, so in other getindex call
@which layer[jmin:jmax, imin:imax] # let's look at getindex(layer::T, i::R, j::R) where {T<:SimpleSDMLayer, R<:UnitRange}
# 2. getindex(layer, jmin:jmax, imin:imax)
T = typeof(layer)
i = jmin:jmax
j = imin:imax
i_min = isempty(i) ? max(i.start-1, 1) : i.start
i_max = isempty(i) ? max(i.stop+2, size(layer, 1)) : i.stop
j_min = isempty(j) ? max(j.start-1, 1) : j.start
j_max = isempty(j) ? max(j.stop+2, size(layer, 2)) : j.stop
i_fix = i_min:i_max
j_fix = j_min:j_max
RT = T <: SimpleSDMResponse ? SimpleSDMResponse : SimpleSDMPredictor
# RT(
    layer.grid[i_fix,j_fix] # ok
    minimum(longitudes(layer)[j_fix])-stride(layer,1) # ok
    maximum(longitudes(layer)[j_fix])+stride(layer,1) # nope
    minimum(latitudes(layer)[i_fix])-stride(layer,2) # nope
    maximum(latitudes(layer)[i_fix])+stride(layer,2) # ok
# )
# Compare allocation for range & LinRange
@time LinRange(layer.bottom+stride(layer, 2), layer.top-stride(layer, 2), size(layer,1))
@time range(layer.bottom+stride(layer, 2), layer.top-stride(layer, 2), length=size(layer,1))
# Let's just switch it to range instead of LinRange
