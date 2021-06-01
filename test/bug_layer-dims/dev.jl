cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.develop(path="$(homedir())/github/SimpleSDMLayers.jl")
using Revise
using SimpleSDMLayers

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

## Step 4: fix for geotiff subset dimensions bug ####

# Load things
layer = SimpleSDMPredictor(WorldClim, BioClim, 1)
layer_subset = layer[coords]
size(layer_subset) # should be 330, 570
left, right, bottom, top = coords
using Test
using ArchGDAL

# Preliminary tests
test = SimpleSDMPredictor(WorldClim, BioClim, 1; coords...)
ref
layer_subset
@test isequal(ref.grid, layer_subset.grid) # true
@test isequal(ref.grid, test.grid) # false
@test size(ref) == size(test) # false
@test stride(ref) == stride(layer_subset) # true
@test stride(ref) == stride(test)

## Investigate function
# SimpleSDMPredictor(...)
resolution = 10.0
file = SimpleSDMLayers._get_raster(WorldClim, BioClim, 1, resolution)
geotiff(SimpleSDMPredictor, file)
geotiff(SimpleSDMPredictor, file; coords...) # wrong dimensions

# geotiff(SimpleSDMPredictor, file; coords...)
import SimpleSDMLayers: _find_span
LT = SimpleSDMPredictor
bandnumber=1
dataset = ArchGDAL.read(file)
# function geotiff(
#     ::Type{LT},
#     file::AbstractString,
#     bandnumber::Integer=1;
#     left = -180.0,
#     right = 180.0,
#     bottom = -90.0,
#     top = 90.0
# ) where {LT<:SimpleSDMLayer}

    # This next block is reading the geotiff file, but also making sure that we
    # clip the file correctly to avoid reading more than we need.
    # This next block is reading the geotiff file, but also making sure that we
    # clip the file correctly to avoid reading more than we need.
    # layer = ArchGDAL.read(file) do dataset

        transform = ArchGDAL.getgeotransform(dataset)
        wkt = ArchGDAL.getproj(dataset)

        # The data we need is pretty much always going to be stored in the first
        # band, but this is not the case for the future WorldClim data.
        band = ArchGDAL.getband(dataset, bandnumber)
        T = ArchGDAL.pixeltype(band)
        
        # The nodata is not always correclty identified, so if it is not found, we assumed it is the smallest value in the band
        nodata = isnothing(ArchGDAL.getnodatavalue(band)) ? convert(T, ArchGDAL.minimum(band)) : convert(T, ArchGDAL.getnodatavalue(band))

        # Get the correct latitudes
        minlon = transform[1]
        maxlat = transform[4]
        maxlon = minlon + size(band,1)*transform[2]
        minlat = maxlat - abs(size(band,2)*transform[6])

        left = isnothing(left) ? minlon : max(left, minlon)
        right = isnothing(right) ? maxlon : min(right, maxlon)
        bottom = isnothing(bottom) ? minlat : max(bottom, minlat)
        top = isnothing(top) ? maxlat : min(top, maxlat)

        lon_stride, lat_stride = transform[2], transform[6]
        
        width = ArchGDAL.width(dataset)
        height = ArchGDAL.height(dataset)

        #global lon_stride, lat_stride
        #global left_pos, right_pos
        #global bottom_pos, top_pos

        lon_stride, left_pos, min_width = _find_span(width, minlon, maxlon, left) # wrong
        _, right_pos, max_width = _find_span(width, minlon, maxlon, right) # ok
        lat_stride, top_pos, max_height = _find_span(height, minlat, maxlat, top) # ok
        _, bottom_pos, min_height = _find_span(height, minlat, maxlat, bottom) # wrong

        max_height, min_height = height .- (min_height, max_height) .+ 1

        # We are now ready to initialize a matrix of the correct type.
        buffer = Matrix{T}(undef, length(min_width:max_width), length(min_height:max_height))
        ArchGDAL.read!(dataset, buffer, bandnumber, min_height:max_height, min_width:max_width)
        buffer = convert(Matrix{Union{Nothing,eltype(buffer)}}, rotl90(buffer))
        replace!(buffer, nodata => nothing)
        LT(buffer, left_pos-0.5lon_stride, right_pos+0.5lon_stride, bottom_pos-0.5lat_stride, top_pos+0.5lat_stride)
    # end

    # return layer

# end

# _find_span
_find_span(width, minlon, maxlon, left) # wrong
_find_span(width, minlon, maxlon, right) # ok
_find_span(height, minlat, maxlat, top) # ok
_find_span(height, minlat, maxlat, bottom) # wrong
n = width
m = minlon
M = maxlon
pos = left
# function _find_span(n, m, M, pos)
    pos > M && return nothing
    pos < m && return nothing
    # stride = (M - m) / n # cannot assign variable to function 🤦
    stride1 = (M - m) / n
    centers = (m + 0.5stride1):stride1:(M-0.5stride1)
    span_pos = last(findmin(abs.(pos .- centers)))
    return (stride, centers[span_pos], span_pos)
# end
# Stupid findmin thing 😡
last(findmin(abs.(pos .- centers)))
pos_diff = abs.(pos .- centers)
last(findmin(pos_diff))
filter(x -> x <= stride1, pos_diff)
findall(x -> x <= stride1, pos_diff)
pos_equal = isapprox.(pos_diff, 0.5stride1)
# Test time
@time _find_span(width, minlon, maxlon, left) # 3 alloc
# With fix
@time _find_span(width, minlon, maxlon, left) # error
@time _find_span(width, minlon, maxlon, left, :none) # error
@time _find_span(width, minlon, maxlon, left, :left) # 10 alloc, ok
@time _find_span(width, minlon, maxlon, right, :right) # ok
@time _find_span(height, minlat, maxlat, top, :top) # ok
@time _find_span(height, minlat, maxlat, bottom, :bottom) # ok
# Test fix
test = SimpleSDMPredictor(WorldClim, BioClim, 1; coords...)
@test size(test) == size(ref)
@test test.grid == ref.grid
@test stride(test) == stride(ref)