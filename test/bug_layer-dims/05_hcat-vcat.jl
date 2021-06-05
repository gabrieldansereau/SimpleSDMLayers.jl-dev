cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.develop(path="$(homedir())/github/SimpleSDMLayers.jl")
using Revise
using SimpleSDMLayers
using Test

## Step 5: investigate error with hcat/vcat
# Checkout previous commit
testpath = "../../../SimpleSDMLayers.jl/test/"
include(joinpath(testpath, "runtests.jl"))
include(joinpath(testpath, "overloads.jl"))

# Failing test
l1 = SimpleSDMPredictor(WorldClim, BioClim, 1; left=0.0, right=10.0, bottom=0.0, top=10.0)
l2 = SimpleSDMPredictor(WorldClim, BioClim, 1; left=0.0, right=10.0, bottom=10.0, top=20.0)
l3 = SimpleSDMPredictor(WorldClim, BioClim, 1; left=10.0, right=20.0, bottom=0.0, top=10.0)
l4 = SimpleSDMPredictor(WorldClim, BioClim, 1; left=10.0, right=20.0, bottom=10.0, top=20.0)

ml1 = hcat(l1, l3)
vl1 = vcat(l1, l2)
ml2 = hcat(l2, l4)
vl2 = vcat(l3, l4)

@test all(vcat(ml1, ml2).grid == hcat(vl1, vl2).grid)

# Keep values
pre = (
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4,
    ml1 = ml1,
    vl1 = vl1,
    ml2 = ml2,
    vl2 = vl2,
)
# Checkout latest commit
new = (
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4,
    ml1 = ml1,
    vl1 = vl1,
    ml2 = ml2,
    vl2 = vl2,
)

# Compare
pre
new
pre.l1[left=0.0, right=10.0, bottom=0.0, top=10.0]
# 
coords_l1 = (left=0.0, right=10.0, bottom=0.0, top=10.0)
coords_l2 = (left=0.0, right=10.0, bottom=10.0, top=20.0)
coords_l3 = (left=10.0, right=20.0, bottom=0.0, top=10.0)
coords_l4 = (left=10.0, right=20.0, bottom=10.0, top=20.0)
boundingbox(l::SimpleSDMLayer) = (left = l.left, right = l.right, bottom = l.bottom, top = l.top)
function fixup(pre, new, coords)
    println(pre)
    println(new)
    fix = pre[coords]
    @test boundingbox(new) == boundingbox(fix)
    @test size(new) == size(fix)
    @test new.grid == fix.grid
    return fix
end
fixup(pre.l1, new.l1, coords_l1)
fixup(pre.l2, new.l2, coords_l2)
fixup(pre.l3, new.l3, coords_l3)
fixup(pre.l4, new.l4, coords_l4)

ml1 = hcat(pre.l1, pre.l3)
ml1 = hcat(new.l1, new.l3)

vl1 = vcat(pre.l1, pre.l2) # whut, no latitude span?
vl1 = vcat(new.l1, new.l2)

ml2 = hcat(pre.l2, pre.l4)
ml2 = hcat(new.l2, new.l4)

vl2 = vcat(pre.l3, pre.l4) # ??
vl2 = vcat(new.l3, new.l4) # ??

vcat(pre.ml1, pre.ml2)
vcat(new.ml1, new.ml2)

hcat(pre.vl1, pre.vl2)
hcat(new.vl1, new.vl2)

# First hcat
l1 = new.l1
l2 = new.l3
t1 = hcat(l1, l2)
boundingbox(t1)
latitudes(t1)
longitudes(t1)

# First vcat
l1 = new.l1
l2 = new.l2
t1 = vcat(l1, l2)
boundingbox(t1)
latitudes(t1)
longitudes(t1)
extrema(latitudes(t1))
Tuple(latitudes(t1)[[1, end]])

coords_ltot = (left = 0.0, right = 20.0, bottom = 0.0, top = 20.0)
ltot = SimpleSDMPredictor(WorldClim, BioClim, 1; coords_ltot...)
ltot.grid
ltot[0.0, 0.0]
ltot[0.0, 20.0]
ltot[20.0, 0.0]
ltot[20.0, 20.0]


# new.l1[0.0, 0.0]
# new.l2[0.0, 20.0]
# new.l3[20.0, 0.0]
# new.l4[20.0, 20.0]