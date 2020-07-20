import Pkg;
Pkg.activate(".")
using SimpleSDMLayers
using DataFrames
using CSV

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
df = CSV.read("data/test-data.csv")
latitude = :latitude
longitude = :longitude
temperature = worldclim(1)
wc_vars = worldclim(1:19)
layer = temperature
layers = wc_vars

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