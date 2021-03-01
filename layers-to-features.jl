import Pkg
Pkg.activate(".")

using Revise
using SimpleSDMLayers
using Plots
using GBIF
using DataFrames

# Get world temperature data
temperature = worldclim(1)

kingfisher = GBIF.taxon("Megaceryle alcyon", strict=true)
kf_occurrences = occurrences(kingfisher)
# Get at least 200 occurrences
while length(kf_occurrences) < 200
    occurrences!(kf_occurrences)
    @info "$(length(kf_occurrences)) occurrences"
end
kf_occurrences

# Clip layer to occurrences
temperature_clip = clip(temperature, kf_occurrences)

presabs1 = mask(temperature_clip, kf_occurrences, Bool)

# Plot occurrences
contour(temperature_clip, fill = true)
scatter!(longitudes(kf_occurrences), latitudes(kf_occurrences))

## Test new mask feature ####

presabs = mask(temperature_clip, kf_occurrences, Bool)
abund = mask(temperature_clip, kf_occurrences, Float64)

presabs.grid
abund.grid

plot(presabs)
plot(abund, c = :BuPu)

presabs2 = convert(Float64, presabs)
presabs.grid
presabs2.grid
abund.grid

tmp = copy(abund)
replace!(tmp.grid, 0 => nothing)
tmp.grid |> unique
abund.grid |> unique
plot(convert(Float64, abund), c = :lightgrey)
plot!(convert(Float64, tmp))

tmp = copy(presabs)
replace!(tmp.grid, 0 => nothing)
tmp.grid |> unique
presabs.grid |> unique
plot(convert(Float64, presabs), c = :lightgrey)
plot!(convert(Float64, tmp))

tmp = copy(presabs)
tmp.grid
tmp = convert(Int, tmp)
tmp.grid
replace!(tmp.grid, 0 => nothing)
tmp.grid |> unique
presabs.grid |> unique
plot(convert(Float64, presabs), c = :lightgrey)
plot!(convert(Float64, tmp))

# Test removezeros argument
presabs1 = mask(temperature_clip, kf_occurrences, Bool)
presabs2 = mask(temperature_clip, kf_occurrences, Bool, removezeros = true)
presabs1.grid |> unique
presabs2.grid |> unique

## Others

grid = [true false false]
grid = [true false NaN]
grid = [true false nothing]
grid = [true false missing]
grid = [1 0 nothing]
grid = [1 0 NaN]

heatmap(grid)

layer = temperature_clip
layer = presabs
lg = copy(layer.grid)

eltype(layer) <: Number

lg = eltype(lg) <: Number ? lg : convert(Matrix{Union{eltype(lg), Float64}}, lg)
replace!(lg, nothing => NaN)
lg = convert(Matrix{Float64}, lg)
longitudes(layer), latitudes(layer), lg


# Sliding window
buffered = slidingwindow(presabs, maximum, 50.0)
plot(buffered)

## Test new replace/replace!

presabs
presabs.grid

tmp = copy(presabs)
replace!(tmp, 0 => nothing)
tmp.grid
tmp.grid |> unique
plot(convert(Float64, tmp), c = :viridis)

## Test DataFrames mask

kf_df = DataFrame(kf_occurrences)
select!(kf_df, [:key, :latitude, :longitude])

# Basic mask

layer = similar(temperature_clip, Bool)
layer = similar(temperature_clip, Float32)

presabs = mask(temperature_clip, kf_df, Bool)
abund   = mask(temperature_clip, kf_df, Float32)

plot(convert(Float32, presabs), c = :viridis)
plot(abund, c = :viridis)

# Latitude/longitude colnames

latitude = :latitude 
longitude = :longitude

kf_df2 = rename(kf_df, :longitude => :lon, :latitude => :lat)
mask(temperature_clip, kf_df2, Bool) # should fail
mask(temperature_clip, kf_df2, Float32) # should fail

mask(temperature_clip, kf_df2, Bool; latitude = :lat, longitude = :lon) # should work
mask(temperature_clip, kf_df2, Float32; latitude = :lat, longitude = :lon) # should work

mask(temperature_clip, kf_df, Bool) # should still work
mask(temperature_clip, kf_df, Float32) # should still work

# Allocations

kf_df
kf_df2 = copy(kf_df)

kf_df2.key[1] = 1
kf_df2
kf_df2[!, :key] .= 1
kf_df2
kf_df

kf_df2 = copy(kf_df)
tmp = mask(temperature_clip, kf_df2, Float32)
tmp.grid
isequal(kf_df, kf_df2)

lat = kf_df[1, :latitude]
lon = kf_df[1, :longitude]
temperature_clip[lon, lat]

tmp = SimpleSDMResponse(kf_df, :key, temperature_clip)
tmp.grid

tempdf = DataFrame(temperature_clip)
templayer = SimpleSDMPredictor(tempdf, :values, temperature_clip)

isequal(templayer, temperature_clip)
isequal(templayer.grid, temperature_clip.grid)

df = tempdf
lats = df[!, latitude]
lons = df[!, longitude]
uniquelats = unique(lats)
uniquelons = unique(lons)
grid = Array{Any}(nothing, size(layer))
lats_idx = [SimpleSDMLayers._match_latitude(layer, lat) for lat in lats]
lons_idx = [SimpleSDMLayers._match_longitude(layer, lon) for lon in lons]
DataFrame(lats = lats_idx, lons = lons_idx)
unique(DataFrame, [:lats, :lons])

## Filter unique sites
tempdf = DataFrame(temperature_clip)
filter!(x -> !isnothing(x.values), tempdf)

# Benchmarks
smallmultidf = vcat(copy(kf_df), copy(kf_df))

mediummultidf = vcat(copy(kf_df), copy(kf_df), copy(kf_df), copy(kf_df),
                     copy(kf_df), copy(kf_df), copy(kf_df), copy(kf_df))

midmediummultidf = vcat(copy(mediummultidf), copy(mediummultidf), copy(mediummultidf),
                        copy(mediummultidf), copy(mediummultidf), copy(mediummultidf))

bigmultidf = vcat(copy(tempdf), copy(tempdf))

hugemultidf = vcat(copy(bigmultidf), copy(bigmultidf), copy(bigmultidf), copy(bigmultidf), 
                   copy(bigmultidf), copy(bigmultidf), copy(bigmultidf), copy(bigmultidf))

# 400 obs (2 x 200)
@time mask(temperature_clip, smallmultidf, Bool) # without unique filtering
# 0.029148 seconds (35.14 k allocations: 5.966 MiB)
# 0.008258 seconds (22.03 k allocations: 5.249 MiB)
@time mask(temperature_clip, smallmultidf, Bool) # with unique filtering
# 0.120525 seconds (75.84 k allocations: 6.536 MiB)
# 0.003245 seconds (12.56 k allocations: 3.843 MiB)

# 1600 obs (8 x 200)
@time mask(temperature_clip, mediummultidf, Bool) # without unique filtering
# 0.038381 seconds (101.17 k allocations: 14.591 MiB)
# 0.030026 seconds (88.04 k allocations: 13.873 MiB, 49.52% gc time)
@time mask(temperature_clip, mediummultidf, Bool) # with unique filtering
# 0.112739 seconds (79.48 k allocations: 6.614 MiB)
# 0.003388 seconds (16.16 k allocations: 3.921 MiB)

# 10 000 obs
@time mask(temperature_clip, midmediummultidf, Bool) # with unique filtering
# 0.118171 seconds (541.20 k allocations: 72.088 MiB, 12.66% gc time)
# 0.092093 seconds (528.03 k allocations: 71.368 MiB, 15.12% gc time)
@time mask(temperature_clip, midmediummultidf, Bool) # with unique filtering
# 0.126579 seconds (103.51 k allocations: 7.160 MiB)
# 0.004906 seconds (40.16 k allocations: 4.466 MiB)

# 100 000 obs
@time mask(temperature_clip, bigmultidf, Bool)
# 0.928445 seconds (7.37 M allocations: 993.580 MiB, 11.33% gc time)
# 0.914822 seconds (7.37 M allocations: 993.580 MiB, 12.57% gc time)
@time mask(temperature_clip, bigmultidf, Bool)
# 0.591735 seconds (3.75 M allocations: 505.812 MiB, 9.40% gc time)
# 0.488477 seconds (3.69 M allocations: 503.119 MiB, 10.87% gc time)

# 1 000 000 obs
@time mask(temperature_clip, hugemultidf, Bool)
# 7.427266 seconds (58.99 M allocations: 7.747 GiB, 11.92% gc time)
# 6.763797 seconds (58.98 M allocations: 7.746 GiB, 10.16% gc time)
@time mask(temperature_clip, hugemultidf, Bool)
# 0.610696 seconds (3.75 M allocations: 528.286 MiB, 9.03% gc time)
# 0.464374 seconds (3.69 M allocations: 525.593 MiB, 8.72% gc time)

# Make sure result is the same
res1 = mask(temperature_clip, smallmultidf, Bool)
res2 = mask(temperature_clip, smallmultidf, Bool)
isequal(res1, res2) # nope, but whatever
isequal(res1.grid, res2.grid) # yup, that's what's important