cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.develop(path="$(homedir())/github/SimpleSDMLayers.jl")
using Revise
using SimpleSDMLayers
using Test

## Step 8 : DataFrames no nothing ####

# For issue
using SimpleSDMLayers
using DataFrames

layer = SimpleSDMPredictor(WorldClim, BioClim, 1)

df = DataFrame([layer, layer])
allowmissing!(df)
for col in [:x1, :x2]
    replace!(df[!, col], nothing => missing)
end
dropmissing(df, [:x1, :x2])

# Tests
DataFrame(layer)
DataFrame([layer, layer])

l1 = layer[top = -89.0]
typeof(l1)
eltype(l1)
typeof(l1.grid)
eltype(l1.grid)

df1 = DataFrame(l1)
typeof(df1.values)
eltype(df1.values)

layers = [layer, layer]
layers = [temperature_clip, temperature_clip]