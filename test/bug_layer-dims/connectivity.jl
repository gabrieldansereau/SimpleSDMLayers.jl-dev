cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.develop(path="$(homedir())/github/SimpleSDMLayers.jl")
using Revise
using SimpleSDMLayers
using Test
using ArchGDAL
using Plots
using StatsBase

file = joinpath(dirname(pathof(SimpleSDMLayers)), "..", "data", "connectivity.tiff")

# geotiff(SimpleSDMPredictor, file)

dataset = ArchGDAL.read(file)
band = ArchGDAL.getband(dataset, 1)
values = ArchGDAL.read(dataset, 1)
ArchGDAL.getgeotransform(dataset)

landcover_mat = permutedims(values[:,end:-1:1])
landcover_mat = convert(Array{Union{Float32, Nothing}}, landcover_mat)
replace!(landcover_mat, -9999.0 => nothing)

coords = (left = -75.17734, right = -72.36486, bottom = 45.34523, top = 47.38457)

mp = SimpleSDMPredictor(landcover_mat, coords...)

qfunc = ecdf(convert(Vector{Float64}, filter(!isnothing, mp.grid)))

qmap = broadcast(qfunc, mp)

plot(qmap, frame=:grid, c=:cork, clim=(0,1))

file2 = joinpath(dirname(pathof(SimpleSDMLayers)), "..", "data", "connectivity2.tiff")
geotiff(file2, mp)
mp2 = geotiff(SimpleSDMPredictor, file2)
mp2 == mp

mp2.grid == mp.grid
mp2.left == mp.left
mp2.right == mp.right
mp2.bottom == mp.bottom
mp2.top == mp.top

# geotiff
dataset = ArchGDAL.read(file2)
band = ArchGDAL.getband(dataset, 1)
values = ArchGDAL.read(dataset, 1)
ArchGDAL.getgeotransform(dataset)

(mp.top - mp.bottom)/size(mp,1)
(mp.right - mp.left)/size(mp,2)

(mp.top - mp.bottom)/2stride(mp, 2)
(mp.right - mp.left)/2stride(mp, 1)

# find_span
layer = copy(mp)
left, right, bottom, top = (-180.0, 180.0, -90.0, 90.0)
bandnumber = 1
dataset = ArchGDAL.read(file2)
import SimpleSDMLayers: _find_span

n, m, M, pos, side = width, minlon, maxlon, right, :right

stride1 = (M - m) / n
centers = (m + 0.5stride1):stride1:(M-0.5stride1)
length(stride1) # length is wrong, should be 1206 (n)

# Whatever, I don't even know if the bounding coordinates are right in the first place...
# Overwriting the tif file
geotiff(file, mp)
