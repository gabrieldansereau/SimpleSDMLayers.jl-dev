import Pkg;
Pkg.activate(".")
using Revise
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