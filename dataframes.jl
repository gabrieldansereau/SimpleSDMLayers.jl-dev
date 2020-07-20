import Pkg;
Pkg.activate(".")
using SimpleSDMLayers
using DataFrames
using CSV

## Prepare test data
#=
df = CSV.read("../betadiversity-hotspots/data/raw/ebd_warblers_head.csv")
newnames = names(df) .|>
    string .|>
    titlecase .|>
    lowercasefirst .|>
    x -> replace(x, " " => "") .|>
    Symbol
rename!(df, newnames)
rename!(df, :scientificName => :species)
df.species = replace.(df.species, " " => "_")
show(df, allcols = true)
newdf = select(df, [:species, :latitude, :longitude])
CSV.write("data/test-data.csv", newdf)
=#

## Test
df = CSV.read("data/test-data.csv")
latitude = :latitude
longitude = :longitude
temperature = worldclim(1)
layer = temperature

# Getindex
temperature[df]

# Setindex
tmp = copy(convert(SimpleSDMResponse, temperature))
tmp[df]
tmp[df] = Float32.(collect(1:nrow(df))) # let's put it aside for now