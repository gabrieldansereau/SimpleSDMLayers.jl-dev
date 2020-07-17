import Pkg
Pkg.activate(".")
using Revise
using SimpleSDMLayers # dev version

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

## Max-min
min(temperature, precipitation)
max(temperature, precipitation)
reduce(max, wc_vars)
wc_vars_resp = convert.(SimpleSDMResponse, wc_vars)
reduce(max, wc_vars_resp)