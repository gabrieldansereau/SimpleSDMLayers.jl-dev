cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.develop(path="$(homedir())/github/SimpleSDMLayers.jl")
using Revise
using SimpleSDMLayers
using Test

## Step 7: No replace! in geotiff ####
import SimpleSDMLayers._prepare_layer_for_burnin
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
layer = SimpleSDMPredictor(WorldClim, BioClim, 1)

array = layer.grid
replace!(array, nothing => NaN)
layer.grid # NaN

array = replace(layer.grid, nothing => NaN)
layer.grid # NaN