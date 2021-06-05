cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.develop(path="$(homedir())/github/SimpleSDMLayers.jl")
using Revise
using SimpleSDMLayers
using Test

# Import QGIS-cropped layer as example
ref = geotiff(SimpleSDMPredictor, "../../rasters/wc1_cropped.tif")

## Step 3: Update match_longitude for range not LinRange ####

# Load things
layer = SimpleSDMPredictor(WorldClim, BioClim, 1)
layer_subset = layer[coords]
size(layer_subset) # should be 330, 570
left, right, bottom, top = coords
import SimpleSDMLayers: _match_latitude, _match_longitude

# Preliminary tests
lon = left
@time _match_longitude(layer, left; side=:left) # 15k alloc
ldiff = abs.(lon .- longitudes(layer))
@time any(x -> isapprox(x, stride(layer, 1)), ldiff)
@time any(x -> isequal(x, stride(layer, 1)), ldiff) # not much difference
@time findlast(x -> isapprox(x, stride(layer, 1)), ldiff)
@time findlast(x -> isequal(x, stride(layer, 1)), ldiff) # not much difference either
@time findfirst(x -> isequal(x, stride(layer, 1)), ldiff)
@time findall(x -> isequal(x, stride(layer, 1)), ldiff)
@time filter(x -> x == stride(layer, 1), ldiff)
@time lon in latitudes(layer)
@time lon in ldiff
@time any(isequal.(ldiff, stride(layer, 1)))
@time any(isapprox.(ldiff, stride(layer, 1)))
@time isequal.(ldiff, stride(layer,1))
# With isequal.
@time _match_longitude(layer, left; side=:left) # 13k
@time _match_longitude(layer, left; side=:none) # 29 alloc (no change)
# with lequal
@time _match_longitude(layer, left; side=:left) # 39 alloc
@time _match_longitude(layer, left; side=:none) # 40 alloc
@time _match_longitude(layer, left; side=:right) # 39 alloc
@time _match_longitude(layer, left+0.001; side=:right) # 41 alloc
# with lapprox
@time _match_longitude(layer, left; side=:left) # 39 alloc
@time _match_longitude(layer, left; side=:none) # 40 alloc
@time _match_longitude(layer, left; side=:right) # 39 alloc
@time _match_longitude(layer, left+0.001; side=:right) # 41 alloc
# Check with latitude
@time _match_latitude(layer, bottom; side=:none) # 40 alloc
@time _match_latitude(layer, bottom; side=:bottom) # 39 alloc
@time _match_latitude(layer, top; side=:top) # 39 alloc
@time _match_latitude(layer, top+0.001; side=:top) # 41 alloc
# Check getindex result
@time layer[coords] # 278 alloc!
# compared to 47.69 k with previous fix, and 222 alloc before previous fix