cd(@__DIR__); import Pkg; Pkg.activate(".")
using SimpleSDMLayers

coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)

## Writing bug as noticed
# Read, write, re-read
ref = SimpleSDMPredictor(WorldClim, BioClim, 1; coords...)
geotiff("test.tif", ref)
test = geotiff(SimpleSDMPredictor, "test.tif")

isequal(size(ref), size(test)) # not equal, 331 x 571 vs 331 x 570

isequal(ref.left, test.left) # equal
isequal(ref.bottom, test.bottom) # equal
isequal(ref.top, test.top) # equal
isequal(ref.right, test.right) # not equal

ref.right # -50.0
test.right # -50.16666666665

# Check the raster file itself
using ArchGDAL
d = ArchGDAL.read("test.tif") # correct raster size, 331 x 571
gt = ArchGDAL.getgeotransform(d)
gt[1] + gt[2]*ArchGDAL.width(d) # correct right bound, ≈ 50.0

# wcpath = joinpath(ENV["SDMLAYERS_PATH"], "WorldClim", "BioClim", "10", "wc2.1_10m_bio_1.tif")
# geotiff(SimpleSDMPredictor, wcpath; coords...)

## Tests in QGIS
# Loaded the world scale temperature data in QGIS
# Cropped it to my coordinates extent
# Dimensions are 330 x 570, not 331 x 571 as I've had for a while
d = ArchGDAL.read("../../rasters/wc1_cropped.tif")
gt = ArchGDAL.getgeotransform(d)
gt[1] + gt[2]*ArchGDAL.width(d)

nlat = (-50.0 - -145.0)/0.1666666666666666574
nlon = (75.0 - 20.0)/0.1666666666666666574

test2 = geotiff(SimpleSDMPredictor, "../../rasters/wc1_cropped.tif")
geotiff("test2.tif", test2)
geotiff(SimpleSDMPredictor, "test2.tif")

# No data value attempt
d = ArchGDAL.read("../../rasters/wc1_cropped.tif")
d1 = ArchGDAL.getband(d, 1)
d1_nodata = ArchGDAL.getnodatavalue(d1)
geotiff("test2.tif", test2; nodata=Float32(-3.4e+38))
geotiff("test2.tif", test2; nodata=Float32(d1_nodata))
geotiff(SimpleSDMPredictor, "test2.tif")
rm("test2.tif")