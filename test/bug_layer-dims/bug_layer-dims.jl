#### Layer dimensions bug

cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)

## Bug as noticed in v0.6.0
cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.add(name = "SimpleSDMLayers", version = "v0.6.0")
using SimpleSDMLayers

t1 = SimpleSDMPredictor(WorldClim, BioClim, 1; coords...) # 331 x 571, as in earlier versions
t2 = SimpleSDMPredictor(WorldClim, BioClim, 1)[coords] # 332 x 571, one extra line!!
isequal(t1, t2) # not equal

## Verification in v0.4.8
# Version used in betadiversity-hotspots before jump to v0.6.0
cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.add(name = "SimpleSDMLayers", version = "v0.4.8")
using SimpleSDMLayers

t1 = worldclim(1; coords...).grid # 331 x 571
t2 = worldclim(1)[coords].grid # 332 x 571!!
isequal(t1, t2) # Nope

## Verification in v0.3.6
# Version used in betadiversity-hotspots before jump to v0.4.8
cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.add(name = "SimpleSDMLayers", version = "v0.3.6")
using SimpleSDMLayers

t1 = worldclim(1; coords...) # 331 x 571
t2 = worldclim(1)[coords] # 331 x 571
isequal(size(t1), size(t2))

## Verification in v0.4.0
cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.add(name = "SimpleSDMLayers", version = "v0.4.0")
using SimpleSDMLayers

t1 = worldclim(1; coords...)
t2 = worldclim(1)[coords]
size(t1) # 331 x 571
size(t2) # 331 x 571
isequal(t1, t2)
# Still works...

## Verification in v0.4.4
# Trying a middle version to find out
cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.add(name = "SimpleSDMLayers", version = "v0.4.4")
using SimpleSDMLayers

t1 = worldclim(1; coords...)
t2 = worldclim(1)[coords]
isequal(size(t1), size(t2))
# Doesn't work

## Verification in v0.4.2
cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.add(name = "SimpleSDMLayers", version = "v0.4.2")
using SimpleSDMLayers

t1 = worldclim(1; coords...)
t2 = worldclim(1)[coords]
isequal(size(t1), size(t2)) # Equal!

isequal(t1, t2) # Nope
# Doesn't work

isequal(t1.grid, t2.grid) # grids are equal
lims1 = t1.left, t1.right, t1.bottom, t1.top
lims2 = t2.left, t2.right, t2.bottom, t2.top
isequal(lims1, lims2) # equal too
# Grids are equal, limits are equal too...
# What's going on??

## Verification in v0.4.3
cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.add(name = "SimpleSDMLayers", version = "v0.4.3")
using SimpleSDMLayers

t1 = worldclim(1; coords...)
t2 = worldclim(1)[coords]
isequal(size(t1), size(t2)) # nope
# So dimensions problem occurred between v0.4.2 and v0.4.3
# But what about equality problem?

## Re-Verification in v0.3.6
cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.add(name = "SimpleSDMLayers", version = "v0.3.6")
using SimpleSDMLayers

t1 = worldclim(1; coords...)
t2 = worldclim(1)[coords]
isequal(size(t1), size(t2))
isequal(t1, t2)
# Ok some equality problem was still around then
# Let's forget about it for now and focus on the dimensions problem
# So back to 0.4.2 and 0.4.3


# Back on v0.4.3
cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.add(name = "SimpleSDMLayers", version = "v0.4.3")
using SimpleSDMLayers

t1 = worldclim(1; coords...)
t2 = worldclim(1)[coords]
isequal(size(t1), size(t2))

t1.grid
t2.grid

t1.grid[end, end]
t2.grid[end, end]

t1.grid[331, 571]
t1[coords.right, coords.top]
t2[coords.right, coords.top]

layer = worldclim(1)
lat_stride = stride(layer, 2)

# _match_longitude
last(findmin(abs.(lon .- longitudes(layer))))

i = SimpleSDMLayers._match_latitude(layer, coords.top)
j = SimpleSDMLayers._match_longitude(layer, coords.right)
layer.grid[i, j]
layer.grid[i-1, j]

last(findmin(abs.(coords.top .- latitudes(layer))))
coords.top .- latitudes(layer)
collect(coords.top .- latitudes(layer))[i]
collect(coords.top .- latitudes(layer))[i-1]
coords.top .- latitudes(layer) |>
    x -> filter(y -> -0.5 < y < 0.5, x)
coords.top .- 0.5*lat_stride .- latitudes(layer) |>
    x -> filter(y -> -0.5 < y < 0.5, x)

i2 = last(findmin(abs.(coords.top .-0.5*lat_stride .- latitudes(layer))))
layer.grid[i2, j]
t1.grid[end, end]

last(findmin(abs.(coords.bottom .- 0.5*lat_stride .- latitudes(layer))))
last(findmin(abs.(coords.bottom .- latitudes(layer))))
coords.bottom .- latitudes(layer) |>
    x -> filter(y -> -0.5 < y < 0.5, x)
coords.bottom .- 0.5*lat_stride .- latitudes(layer) |>
    x -> filter(y -> -0.5 < y < 0.5, x)

last(findmin(abs.(coords.right .- 0.5*lat_stride .- longitudes(layer))))
last(findmin(abs.(coords.right .- longitudes(layer))))
last(findmin(abs.(coords.right .- longitudes(layer))))

## Dev some changes
cd(@__DIR__); import Pkg; Pkg.activate(".")
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
Pkg.develop(path="$(homedir())/github/SimpleSDMLayers.jl")
using Revise
using SimpleSDMLayers

t1 = SimpleSDMPredictor(WorldClim, BioClim, 1; coords...) # 331 x 571, as in earlier versions
t2 = SimpleSDMPredictor(WorldClim, BioClim, 1)[coords] # 332 x 574, one extra line!!
isequal(t1, t2) # not equal

using Test
S = SimpleSDMResponse(rand(Bool, 4, 4), 0.3, 0.7, 0.5, 0.9)
for (i,l) in enumerate(latitudes(S))
    @test SimpleSDMLayers._match_latitude(S, l) == i
end
[SimpleSDMLayers._match_latitude(S, l) for l in latitudes(S)]

for (i,l) in enumerate(longitudes(S))
    @test SimpleSDMLayers._match_longitude(S, l) == i
end
[SimpleSDMLayers._match_latitude(S, l) for l in latitudes(S)]