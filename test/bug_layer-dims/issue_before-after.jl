cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.develop(path="$(homedir())/github/SimpleSDMLayers.jl")
using Revise
using SimpleSDMLayers
using Test

## Step 7: Final before & after

temp = SimpleSDMPredictor(WorldClim, BioClim, 1);
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0);

## Before
# getindex
l1 = temp[coords];
size(l1)
(l1.bottom, l1.top)
(l1.left, l1.right)
stride(l1)

# geotiff reading
l2 = SimpleSDMPredictor(WorldClim, BioClim, 1; coords...);
l2.grid == l1.grid
size(l2)
(l2.bottom, l2.top)
(l2.left, l2.right)
stride(l2)

# geotiff writing
l2_backup = copy(l2);
tempfile = tempname();
geotiff(tempfile, l2);
l3 = replace(geotiff(SimpleSDMPredictor, tempfile), NaN => nothing);
l2_backup.grid == l2.grid
l2_backup.grid == l3.grid

# Inverse bounds
SimpleSDMPredictor(temp.grid, 180.0, -180.0, 90.0, -90.0)

# hcat/vcat
l4 = SimpleSDMPredictor(WorldClim, BioClim, 1; left=0.0, right=10.0, bottom=0.0, top=10.0);
l5 = SimpleSDMPredictor(WorldClim, BioClim, 1; left=0.0, right=10.0, bottom=10.0, top=20.0);
l6 = SimpleSDMPredictor(WorldClim, BioClim, 1; left=0.0, right=10.0, bottom=20.0, top=30.0);
vcat(l4, l5)
vcat(l4, l6)