## Run tests in proper folder

# Go to test repo
cd("$(homedir())/github/SimpleSDMLayers.jl/test")
import Pkg; Pkg.activate(".")

# Change package version
Pkg.add(name = "DataFrames", version = "v1.1")
Pkg.add(name = "DataFrames", version = "v1.0")
Pkg.add(name = "DataFrames", version = "v0.22")
Pkg.add(name = "DataFrames", version = "v0.22.0")
Pkg.add(name = "DataFrames", version = "v0.21")
# Run tests
include(abspath("runtests.jl"))

## Fix problem with DataFrames v1.0
cd("$(homedir())/github/SimpleSDMLayers.jl/test")
import Pkg; Pkg.activate(".")

# Change package version
Pkg.add(name = "DataFrames", version = "v1.0")
Pkg.add(name = "DataFrames", version = "v0.22")
Pkg.add(name = "DataFrames", version = "v0.22.0")
Pkg.add(name = "DataFrames", version = "v0.21")

using Revise
using SimpleSDMLayers
using DataFrames

temperature = SimpleSDMPredictor(WorldClim, BioClim, 1)
df = DataFrame(latitude = [0.0, 1.0], longitude = [30.0, 31.0], values = [42.0, 15.0])
temperature_clip = clip(temperature, df)
typeof(DataFrame([temperature_clip, temperature_clip])) == DataFrame

# Single layer conversion
DataFrame(temperature_clip)

# Multiple layer conversion
layers = [temperature_clip, temperature_clip]
DataFrame(layers)
# DataFrame(layers, :auto)

## Test GBIF compatibility
cd("$(homedir())/github/SimpleSDMLayers.jl/test")
import Pkg; Pkg.activate(".")

Pkg.add(name = "GBIF", version = "v0.2")
Pkg.add(name = "GBIF", version = "v0.3")

include(abspath("gbif.jl"))
include(abspath("runtests.jl"))

## Test Plots compatibility
cd("$(homedir())/github/SimpleSDMLayers.jl/test")
import Pkg; Pkg.activate(".")

Pkg.add(name = "Plots", version = "v1.14.0")
Pkg.add(name = "Plots", version = "v1.9.0")
Pkg.add(name = "Plots", version = "v1.0.0") # fails because of color :terrain
Pkg.add(name = "Plots", version = "v1.5.0")
Pkg.add(name = "Plots", version = "v1.2.0")
Pkg.add(name = "Plots", version = "v1.1.0")
Pkg.add(name = "Plots", version = "v1.0")
Pkg.add(name = "Plots", version = "v0")

include(abspath("plots.jl"))
include(abspath("runtests.jl"))

temperature, precipitation = SimpleSDMPredictor(WorldClim, BioClim, [1,12])
plot(precipitation)
histogram(precipitation, leg=false)

## Test main compat entries
cd("$(homedir())/github/SimpleSDMLayers.jl/test")
import Pkg; Pkg.activate(".")

Pkg.add(name = "ArchGDAL", version = "v0.4") # fails in overloads.jl, something with size
Pkg.add(name = "ArchGDAL", version = "v0.5") # fails in chelsa.jl, something with Int16 in RCP26
Pkg.add(name = "ArchGDAL", version = "v0.6")
include(abspath("runtests.jl"))

Pkg.add(name = "Downloads", version = "v1.4")
Pkg.add(name = "Downloads", version = "v1.3") # unsatistiable requirements
include(abspath("runtests.jl"))

Pkg.add(name = "RecipesBase", version = "v0.7")
Pkg.add(name = "RecipesBase", version = "v0.8") # fails, but v0.7 works?? # nope, worked wthis time
Pkg.add(name = "RecipesBase", version = "v1.0.0")
include(abspath("runtests.jl"))

Pkg.add(name = "Requires", version = "v1.0")
include(abspath("runtests.jl"))

Pkg.add(name = "ZipFile", version = "v0.8")
Pkg.add(name = "ZipFile", version = "v0.9")
include(abspath("runtests.jl"))
