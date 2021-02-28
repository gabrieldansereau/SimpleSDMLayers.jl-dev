import Pkg
Pkg.activate(".")

using Revise
using SimpleSDMLayers
using Plots
using GBIF

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