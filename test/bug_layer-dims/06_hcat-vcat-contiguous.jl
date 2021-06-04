cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.develop(path="$(homedir())/github/SimpleSDMLayers.jl")
using Revise
using SimpleSDMLayers
using Test

## Step 6: Contiguous hcat/vcat ####
using Test

# Base layer
l1 = SimpleSDMPredictor(WorldClim, BioClim, 1; left=0.0, right=10.0, bottom=0.0, top=10.0)
# Non contiguous layers
l2 = SimpleSDMPredictor(WorldClim, BioClim, 1; left=0.0, right=10.0, bottom=20.0, top=30.0)
l3 = SimpleSDMPredictor(WorldClim, BioClim, 1; left=20.0, right=30.0, bottom=0.0, top=10.0)

vcat(l2, l1)
hcat(l1, l3)

l4 = SimpleSDMPredictor(WorldClim, BioClim, 1; left=0.0, right=10.0, bottom=10.0, top=20.0)
vcat(l4, l1)
vcat(l1, l4)

## Step 7: No replace! in geotiff ####
import SimpleSDMLayers._prepare_layer_for_burnin
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
layer = SimpleSDMPredictor(WorldClim, BioClim, 1)

array = layer.grid
replace!(array, nothing => NaN)
layer.grid # NaN

array = replace(layer.grid, nothing => NaN)
layer.grid # NaN