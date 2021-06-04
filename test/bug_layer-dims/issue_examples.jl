## Issue examples ####

# Dimensions should be 330 x 570
cd("assets/")

wcpath = joinpath(ENV["SDMLAYERS_PATH"], "WorldClim", "BioClim", "10", "wc2.1_10m_bio_1.tif")
tmpfile = tempname()

query = `gdalwarp -te -145.0 20.0 -50.0 75.0 $(wcpath) $(tmpfile)`
run(query)

using ArchGDAL
d = ArchGDAL.read(tmpfile)
ArchGDAL.height(d) # 330
ArchGDAL.width(d) # 570

# hcat/vcat are wrong
l1 = SimpleSDMPredictor(WorldClim, BioClim, 1; left=0.0, right=10.0, bottom=0.0, top=10.0);
l2 = SimpleSDMPredictor(WorldClim, BioClim, 1; left=0.0, right=10.0, bottom=10.0, top=20.0);
l3 = SimpleSDMPredictor(WorldClim, BioClim, 1; left=10.0, right=20.0, bottom=0.0, top=10.0);

vl1 = vcat(l1, l2) # latitudes span makes no sense
ml1 = hcat(l1, l3) # things look fine here, but are they?

latitudes(vl1)
vl1.bottom < vl1.top # false

latitudes(ml1)
ml1.bottom < ml1.top # false

# Inversed bounds
l1 = SimpleSDMPredictor(WorldClim, BioClim, 1)
l2 = SimpleSDMPredictor(copy(l1.grid), 180.0, -180.0, 90.0, -90.0) # looks the same, but bounds really are inversed

l1.grid == l2.grid # same grid
extrema(longitudes(l1)) == extrema(longitudes(l2)) # same longitude extremas (displayed by show method)
extrema(latitudes(l1)) == extrema(latitudes(l2)) # same latitude extremas too

longitudes(l1) == longitudes(l2) # false, inversed
latitudes(l1) == latitudes(l2) # false, inversed

# Geotiff writing
l1 = SimpleSDMPredictor(WorldClim, BioClim, 1)
length(l1)
geotiff(tempname(), l1)
length(l1)