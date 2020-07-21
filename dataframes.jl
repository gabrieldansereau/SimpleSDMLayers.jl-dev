import Pkg;
Pkg.activate(".")
using SimpleSDMLayers
using DataFrames
using CSV
using Plots

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
testdf = CSV.read("data/test-data.csv")
latitude = :latitude
longitude = :longitude
temperature = worldclim(1)
wc_vars = worldclim(1:19)
layer = temperature
layers = wc_vars
col = :values
ty = :SimpleSDMResponse

# Getindex
temperature[df]

# Setindex
tmp = copy(convert(SimpleSDMResponse, temperature))
tmp[df]
# tmp[df] = Float32.(collect(1:nrow(df))) # let's put it aside for now

# Inconsistent use of missing/nothing/NaN in GBIF code
using GBIF
kingfisher = GBIF.taxon("Megaceryle alcyon", strict=true)
kf_occurrences = occurrences(kingfisher)
for i in 1:9
  occurrences!(kf_occurrences)
end
tmp = [k.latitude for k in kf_occurrences]
sort(tmp) # none missing
tmp = layer[kf_occurrences]
sort(tmp) # none missing
temperature_clip = clip(temperature, kf_occurrences)
# whatever...

# Clip
temperature_clip = clip(temperature, df)
temperature_clip.grid
# using Plots
# plot(temperature_clip)

# DataFrame
@time DataFrame(temperature)
@time DataFrame(layers)

# SimpleSDMResponse/Predictor
df = DataFrame(temperature)
@time SimpleSDMResponse(df, :values, layer)
@time SimpleSDMPredictor(df, :values, layer)
minidf = DataFrame(latitude = [1.0, 0.0], longitude = [1.0, 0.0], values = [42, 3])
@time SimpleSDMResponse(minidf, :values, layer)
mediumdf = df[rand(1:nrow(df), 10_000), :]
@time SimpleSDMResponse(mediumdf, :values, layer)
valuesdf = filter(x -> !isnothing(x.values), df)
@time tmp = SimpleSDMResponse(valuesdf, :values, layer)
