cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.develop(path="$(homedir())/github/SimpleSDMLayers.jl")
using Revise
using SimpleSDMLayers
using Test

# Import QGIS-cropped layer as example
ref = geotiff(SimpleSDMPredictor, "../../rasters/wc1_cropped.tif")

## Step 1: Fix getindex ####

layer = SimpleSDMPredictor(WorldClim, BioClim, 1)
layer_subset = layer[coords]
size(layer_subset) # should be 330, 570

left, right, bottom, top = coords

import SimpleSDMLayers: _match_latitude, _match_longitude

# Inside getindex call
imax = _match_longitude(layer, isnothing(right) ? layer.right : right)
imin = _match_longitude(layer, isnothing(left) ? layer.left : left)
jmax = _match_latitude(layer, isnothing(top) ? layer.top : top)
jmin = _match_latitude(layer, isnothing(bottom) ? layer.bottom : bottom)
test1a = layer[jmin:jmax, imin:imax]
test1b = layer[jmin+1:jmax-1, imin+1:imax] # almost same as ref
ref

isequal(ref.grid, test1b.grid) # equal

stride(ref)
stride(layer) # original stride ok
stride(test1a) # stride is wrong...
stride(test1b) # still wrong

longitudes(ref)
longitudes(test1a)
longitudes(test1b)
isequal(collect(longitudes(ref)), collect(longitudes(test1b)))
isequal(collect(longitudes(ref)), collect(longitudes(test1b)))

latitudes(test1)
latitudes(ref)

# Check functions 1 by 1
# getindex(layer; coords...) # should be ok in theory
# _match_longitude # could be the problem
_match_longitude(layer, isnothing(right) ? layer.right : right) # ok
_match_longitude(layer, isnothing(left) ? layer.left : left) # should be 211
lon = left
last(findmin(abs.(lon .- longitudes(layer))))
@time lon .- longitudes(layer)
@time abs.(lon .- longitudes(layer))
@time findmin(abs.(lon .- longitudes(layer)))
@time last(findmin(abs.(lon .- longitudes(layer))))

isapprox(abs.(lon .- longitudes(layer))[[210, 211]]...)
findall(x -> x ≈ stride(layer,1), abs.(lon .- longitudes(layer)))

@time filter(x -> lon-stride(layer,1) <= x <= lon+stride(layer, 1), longitudes(layer))
i1, i2 = findall(x -> lon-stride(layer,1) <= x <= lon+stride(layer, 1), longitudes(layer))
(lon .- longitudes(layer))[[i1, i2]] .|> abs
isapprox.(stride(layer,1), abs.((lon .- longitudes(layer))[[i1, i2]]))

@time any(x -> isapprox(x, left), longitudes(layer) .- stride(layer,1))
findall(x -> isapprox(x, left), longitudes(layer) .- stride(layer,1))
findall(x -> isapprox(x, left), longitudes(layer) .+ stride(layer,1))
findall(x -> isapprox(x, left+0.001), longitudes(layer) .- stride(layer,1))
findall(x -> isapprox(x, right), longitudes(layer) .- stride(layer,1))

findall(x -> isapprox(x, bottom), latitudes(layer) .- stride(layer,1))
findall(x -> isapprox(x, top), latitudes(layer) .+ stride(layer,1))


@time LinRange(layer.bottom+stride(layer, 2), layer.top-stride(layer, 2), size(layer,1))
@time range(layer.bottom+stride(layer, 2), layer.top-stride(layer, 2), length=size(layer,1))

# As previously
@time _match_longitude(layer, isnothing(left) ? layer.left : left) # 25 alloc
# with return l
@time _match_longitude(layer, isnothing(left) ? layer.left : left) # 25 alloc
# with if else
@time _match_longitude(layer, isnothing(left) ? layer.left : left) # 25 alloc
@time _match_longitude(layer, isnothing(left) ? layer.left : left; side=:none) # 25 alloc
@time _match_longitude(layer, isnothing(left) ? layer.left : left; side=:left) # 25 alloc
@time _match_longitude(layer, isnothing(left) ? layer.left : left; side=:youppi) # 25 alloc
# with any statement
ldiff = abs.(lon .- longitudes(layer))
@time any(x -> isapprox(x, stride(layer,1)), ldiff)
@time findlast(x -> isapprox(x, stride(layer,1)), ldiff)
@time _match_longitude(layer, isnothing(left) ? layer.left : left) # 26 alloc
@time _match_longitude(layer, isnothing(left) ? layer.left : left; side=:none) # 26 alloc
@time _match_longitude(layer, isnothing(left) ? layer.left : left; side=:left) # 15.15k alloc
@time _match_longitude(layer, isnothing(right) ? layer.right : right; side=:right) # 10.95k alloc
@time _match_longitude(layer, isnothing(left+0.001) ? layer.left : left+0.001; side=:left) # 15.15k alloc
@time _match_longitude(layer, isnothing(left-0.001) ? layer.left : left-0.001) # 28 alloc
@time _match_longitude(layer, layer.right) # 28 alloc
@time _match_longitude(layer, layer.right; side=:right) # 28 alloc
@time _match_longitude(layer, layer.right; side=:left) # 28 alloc
# Test full fix
_match_longitude(layer, isnothing(right) ? layer.right : right; side=:right)
_match_longitude(layer, isnothing(left) ? layer.left : left; side=:left)
_match_latitude(layer, isnothing(top) ? layer.top : top; side=:top)
_match_latitude(layer, isnothing(bottom) ? layer.bottom : bottom; side=:bottom)
@time layer[coords] # 0.002510 seconds (47.69 k allocations: 1.676 MiB) (now)
@time layer[coords] # 0.001570 seconds (222 allocations: 981.422 KiB) (before)

ref
isequal(layer[coords].grid, ref.grid) # true!
isequal(layer[coords], ref) # false, as bounding coordinates are not exactly the same

tmp = collect(range(layer.left+stride(layer, 1), layer.right-stride(layer, 1), length=size(layer,2)))
145.0 in tmp .+ stride(layer,1)
# longitudes # should be ok in theory
# stride # should be ok in theory
# getindex(layer, irange, jrange) # could be the problem too
i = jmin:jmax # need to be inversed
j = imin:imax
layer.grid[i, j] # wrong dimensions
i_fix = i
j_fix = j
layer.grid[i_fix,j_fix]
minimum(longitudes(layer)[j_fix])-stride(layer,1)
maximum(longitudes(layer)[j_fix])+stride(layer,1)
minimum(latitudes(layer)[i_fix])-stride(layer,2)
maximum(latitudes(layer)[i_fix])+stride(layer,2)