## Some tests
import Pkg; Pkg.activate(".")
using Revise
using SimpleSDMLayers
using DataFrames

# Generate some elements
testvec = [1, 2, 3]
testmat = hcat(testvec, testvec)
testdf = DataFrame(testmat, :auto)
testlayer = SimpleSDMPredictor(testmat)
# Generate some arrays of elements
arrayvec = [testvec for i in 1:15]
arraymat = [testmat for i in 1:15]
arraydf = [testdf for i in 1:15]
arraylayer = [testlayer for i in 1:15]
arraymix = [testvec, testmat, testdf, testlayer]
# And a tuple for fun
(testvec, testmat, testdf, testlayer)

# Show the layer types
testlayer
arraylayer

layer = SimpleSDMPredictor(WorldClim, BioClim, 1)
layers = SimpleSDMPredictor(WorldClim, BioClim, 1:2)

layer
layers


## Learn about custom pretty printing
# https://docs.julialang.org/en/v1/manual/types/#man-custom-pretty-printing
# https://discourse.julialang.org/t/overload-show-for-array-of-custom-types/9589
# https://discourse.julialang.org/t/extending-base-show-for-array-of-types/31289
# https://stackoverflow.com/questions/58962304/how-to-overload-base-show-for-custom-array-types

## Play with show overloads
# Original overload
function Base.show(io::IO, layer::T) where {T <: SimpleSDMLayer}
    itype = eltype(layer)
    otype = T <: SimpleSDMPredictor ? "predictor" : "response"
    print(io, """SDM $(otype) → $(size(layer,1))×$(size(layer,2)) grid with $(length(layer)) $(itype)-valued cells
    Latitudes\t$(extrema(latitudes(layer)))
    Longitudes\t$(extrema(longitudes(layer)))""")
end
testlayer
arraylayer
# Spaces to indent coordinates
function Base.show(io::IO, layer::T) where {T <: SimpleSDMLayer}
    itype = eltype(layer)
    otype = T <: SimpleSDMPredictor ? "predictor" : "response"
    print(io, """SDM $(otype) → $(size(layer,1))×$(size(layer,2)) grid with $(length(layer)) $(itype)-valued cells
    \x20\x20Latitudes\t$(extrema(latitudes(layer)))
    \x20\x20Longitudes\t$(extrema(longitudes(layer)))""")
end
testlayer
arraylayer

# Spaces to indent coordinates, but slightly different print call
function Base.show(io::IO, layer::T) where {T <: SimpleSDMLayer}
    itype = eltype(layer)
    otype = T <: SimpleSDMPredictor ? "predictor" : "response"
    print(io, 
          "SDM $(otype) → $(size(layer,1))×$(size(layer,2)) grid with $(length(layer)) $(itype)-valued cells\n",
          "\x20\x20Latitudes\t$(extrema(latitudes(layer)))\n",
          "\x20\x20Longitudes\t$(extrema(longitudes(layer)))"
          )
end
testlayer
arraylayer
## Different show for single element vs in an Array
# Single
function Base.show(io::IO, ::MIME"text/plain", layer::T) where {T <: SimpleSDMLayer}
    itype = eltype(layer)
    otype = T <: SimpleSDMPredictor ? "predictor" : "response"
    print(io, """SDM $(otype) → $(size(layer,1))×$(size(layer,2)) grid with $(length(layer)) $(itype)-valued cells
    \x20\x20Latitudes\t$(extrema(latitudes(layer)))
    \x20\x20Longitudes\t$(extrema(longitudes(layer)))""")
end
testlayer
# In an Array and probably everywhere else
function Base.show(io::IO, layer::T) where {T <: SimpleSDMLayer}
    itype = eltype(layer)
    otype = T <: SimpleSDMPredictor ? "predictor" : "response"
    print(io, "SDM $(otype) → $(size(layer,1))×$(size(layer,2)) grid with $(length(layer)) $(itype)-valued cells")
end
arraylayer

function Base.show(io::IO, ::MIME"text/plain", layer::T) where {T <: SimpleSDMLayer}
    itype = eltype(layer)
    otype = T <: SimpleSDMPredictor ? "predictor" : "response"
    print(io, """SDM $(otype) → $(size(layer,1))×$(size(layer,2)) grid with $(length(layer)) $(itype)-valued cells
    \x20\x20Latitudes\t$(extrema(latitudes(layer)))
    \x20\x20Longitudes\t$(extrema(longitudes(layer)))""")
end