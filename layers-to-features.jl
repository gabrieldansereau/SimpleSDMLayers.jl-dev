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