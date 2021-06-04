cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.develop(path="$(homedir())/github/SimpleSDMLayers.jl")
using Revise
using SimpleSDMLayers
using Test

## Step 9 : Minor improvements ####

## 1. geotiff(file; coords)
file = "assets/wc2.1_10m_bio_1.tif"
methods(geotiff)
geotiff(SimpleSDMPredictor, file, 1; coords...)

geotiff(SimpleSDMPredictor, file, 1, coords)
geotiff(SimpleSDMPredictor, file, 1, coords)
# Whatever, too many methods to define and too hard to maintain

## 2. ==
l1, l2 = SimpleSDMPredictor(WorldClim, BioClim, 1:2)
l3 = copy(l1)
l4 = similar(l1)
replace!(l4, nothing => NaN)
l5 = SimpleSDMPredictor(replace(l1.grid, nothing => missing), l1)

@test l1 == l1
@test l1 === l1
@test l2 != l1
@test l3 == l1
@which l3 !== l1

@test l4 != l1
@test l4 != l4
@test !isequal(l4, l1)
@test isequal(l4, l4)

@test ismissing(l5 == l1)
@test ismissing(l5 == l5)
@test !isequal(l5, l1)
@test isequal(l5, l5)

# basepath = expanduser("~/julia/julia-1.6.1/share/julia/base/")
# @which isequal(l3, l1) # Base.operators.jl:123
# run(`code $(basepath)/operators.jl`)
# @which l3 == l1 # Base.Base.jl:87
# run(`code $(basepath)/Base.jl`)
# @which isequal(l3.grid, l1.grid) # Base.abstractarray.jl:1962
# @which l3.grid == l1.grid # Base.abstractarray.jl:1991
# run(`code $(basepath)/abstractarray.jl`)

# 3. replace!
l1 = SimpleSDMPredictor(WorldClim, BioClim, 1)
l2 = convert(SimpleSDMResponse, l2)
@which replace!(l1, nothing => NaN)
l1 isa SimpleSDMPredictor
l2 isa SimpleSDMPredictor
replace!(l1, nothing => NaN)
replace!(l2, nothing => NaN)
l2.grid

## 4. layer1[layer2]

l1, l2 = SimpleSDMPredictor(WorldClim, BioClim, 1:2)
l3 = SimpleSDMPredictor(WorldClim, BioClim, 1; left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
stride(l1) == stride(l3) # same stride

l4 = l1[l2] # works
l4 == l1 # but it's exactly the same

l5 = l1[l3] # doesn't work as layers have different sizes, but should work
l5 == l3 # should return true

l3[l1]

## 5. boundingbox(layer)

## 6. hash() method for custom type

l1 = SimpleSDMPredictor(WorldClim, BioClim, 1; left = 0.0, right = 10.0, bottom = 0.0, top = 10.0)
l2 = copy(l1)

l1 == l2
hash(l1) == hash(l2)
hash(l1.grid) == hash(l2.grid)

## 7. WorldClim no data
using ArchGDAL
import SimpleSDMLayers: _find_span
import SimpleSDMLayers: _prepare_layer_for_burnin
coords = (left = 0.0, right = 10.0, bottom = 0.0, top = 10.0)
layer = SimpleSDMPredictor(WorldClim, BioClim, 1; coords...)

wcfile = "assets/wc2.1_10m_bio_1.tif"
tempfile = tempname()
bandnumber = 1
left, right, bottom, top = fill(nothing, 4)

geotiff(tempfile, layer)
geotiff(SimpleSDMPredictor, tempfile)
geotiff(SimpleSDMPredictor, wcfile)

dataset = ArchGDAL.read(tempfile)
band = ArchGDAL.getband(dataset, bandnumber)
ArchGDAL.getnodatavalue(band)

dataset = ArchGDAL.read(wcfile)
band = ArchGDAL.getband(dataset, bandnumber)
nd = ArchGDAL.getnodatavalue(band)
nd2 = Float32(nd)

geotiff(tempfile, layer; nodata = Float32(nd))
geotiff(tempfile, layer; nodata = Float32(nd2))
geotiff(SimpleSDMPredictor, tempfile)

l = SimpleSDMPredictor(WorldClim, BioClim, 1)

f2 = tempname()
geotiff(f2, l; nodata=-3.4f38)
mp2 = geotiff(SimpleSDMResponse, f2)

@test typeof(mp2) <: SimpleSDMResponse
@test size(mp2) == size(l)
@test mp2 == l

## Issue comment
l = SimpleSDMPredictor(WorldClim, BioClim, 1)

# With -9999 as nodata
f1 = tempname()
geotiff(f1, l) # use -9999 as nodata by default
mp1 = geotiff(SimpleSDMResponse, f1)

length(mp1) == length(l) # false, but should be true
mp1 == l # false, same thing
replace(mp1, NaN => nothing) == l # true, but could be avoided

# With custom nodata
f2 = tempname()
geotiff(f2, l; nodata=-3.4f38) # this is the nodata from the WorldClim files
mp2 = geotiff(SimpleSDMPredictor, f2)

length(mp2) == length(l) # false
mp2 == l # false
replace(mp2, NaN => nothing) == l # true
