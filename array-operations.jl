import Pkg
Pkg.activate(".")
using Revise
using SimpleSDMLayers # dev version
using Statistics

temperature = worldclim(1)
precipitation = worldclim(12)
wc_vars = worldclim(1:19)

# Extrema overload
minimum(temperature)
maximum(temperature)
extrema(temperature)

# Broadcast type
broadcast(abs, temperature)
maximum(abs(temperature))

# Max-min
@time min(temperature, precipitation)
@time max(temperature, precipitation)
@time test = reduce(max, wc_vars)
wc_vars_resp = convert.(SimpleSDMResponse, wc_vars)
@time test = reduce(max, wc_vars_resp)

minifloat = SimpleSDMPredictor([-1 0 1 nothing])
miniint = SimpleSDMResponse([-1.0 0.0 1.0 nothing])
mini32 = SimpleSDMPredictor([Float32.([-1.0 0.0 1.0])... nothing])
min(miniint, minifloat)
@time reduce(min, [miniint, minifloat, mini32])

# Mean-std
mean(wc_vars)
std(wc_vars)
# Attempt 1
mininothing = SimpleSDMPredictor([nothing 1.0 nothing nothing])
coordsnothing = map(x -> findall(isnothing, x.grid), [minifloat, miniint, mini32, mininothing])
uniquenothing = unique(coordsnothing)
unique(reduce(vcat, uniquenothing)) |> sort
# Need InvertedIndex...
# Attempt 2
coordsnotnothing = map(x -> findall(!isnothing, x.grid), [minifloat, miniint, mini32, mininothing])
uniquenot = unique(coordsnotnothing)
@time reduce(intersect, uniquenot)
@time intersect(uniquenot...)
