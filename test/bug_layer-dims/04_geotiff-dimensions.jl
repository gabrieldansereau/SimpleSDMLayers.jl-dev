cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.develop(path="$(homedir())/github/SimpleSDMLayers.jl")
using Revise
using SimpleSDMLayers
using Test

# Import QGIS-cropped layer as example
ref = geotiff(SimpleSDMPredictor, "../../rasters/wc1_cropped.tif")

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