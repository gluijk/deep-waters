# Speed comparison raster vs terra packages
# www.overfitting.net
# https://www.overfitting.net/2023/12/simulacion-de-diluvio-universal-con-r.html

library(raster)  # read GeoTIFF, reprojection, crop and resample
library(terra)  # read GeoTIFF, reprojection, crop and resample


# READ raster
deepwaters=raster("gebco_2023_n69.5654_s25.1807_w-25.752_e31.2891.tif")  # read GeoTIFF file
deepwaters
plot(deepwaters)
abline(v=0)  # Greenwich meridian

deepwaters2=rast("gebco_2023_n69.5654_s25.1807_w-25.752_e31.2891.tif")
deepwaters2
plot(deepwaters2)
abline(v=0)  # Greenwich meridian


# REPROJECT raster from Longitude Latitude (+proj=longlat)/WGS84
# to Lambert Conic Conformal (+proj=lcc)/WGS84
start_time=Sys.time()
CRS="+proj=lcc +ellps=WGS84 +lat_1=33 +lat_2=45 +lon_0=0 +units=km"
deepwatersrp=projectRaster(deepwaters, crs=CRS)
end_time=Sys.time()
lapso_raster=end_time-start_time  # Time difference of 4.666385 mins
deepwatersrp

start_time=Sys.time()
CRS="+proj=lcc +ellps=WGS84 +lat_1=33 +lat_2=45 +lon_0=0 +units=km"
deepwatersrp2=project(x=deepwaters2, y=CRS)
end_time=Sys.time()
lapso_terra=end_time-start_time  # Time difference of 1.858117 mins
deepwatersrp2


# PLOT raster
start_time=Sys.time()
for (i in 1:10) {
    plot(deepwatersrp)
    abline(v=0)  # Greenwich meridian
}
end_time=Sys.time()
lapso_raster=end_time-start_time  # Time difference of 1.718697 mins

start_time=Sys.time()
for (i in 1:10) {
    plot(deepwatersrp2)
    abline(v=0)  # Greenwich meridian
}
end_time=Sys.time()
lapso_terra=end_time-start_time  # Time difference of 10.92895 secs


# CROP raster to area of interest (save Canary Islands and The Netherlands)
start_time=Sys.time()
cropdef=as(extent(-1850, 2200, 3450, 6050), 'SpatialPolygons')  # kms
crs(cropdef)=crs(deepwatersrp)
deepwaterscrop=crop(deepwatersrp, cropdef)
end_time=Sys.time()
lapso_raster=end_time-start_time  # Time difference of 4.796154 secs
deepwaterscrop
plot(deepwaterscrop)
abline(v=0)  # Greenwich meridian

start_time=Sys.time()
cropdef2=ext(-1850, 2200, 3450, 6050)
deepwaterscrop2=crop(deepwatersrp2, cropdef2)
end_time=Sys.time()
lapso_terra=end_time-start_time  # Time difference of 2.65933 secs
deepwaterscrop2
plot(deepwaterscrop2)
abline(v=0)  # Greenwich meridian


# RESAMPLE raster to Full HD
DIMY=1080
DIMX=1920

start_time=Sys.time()
deepwatersrs=raster(nrow=DIMY, ncol=DIMX, ext=extent(deepwaterscrop), crs=crs(deepwaterscrop))
deepwatersrs=resample(deepwaterscrop, deepwatersrs, method='bilinear')  # bilinear is faster than ngb!
end_time=Sys.time()
lapso_raster=end_time-start_time  # Time difference of 25.8529 secs
plot(deepwatersrs)
abline(v=0)  # Greenwich meridian

# https://www.pmassicotte.com/posts/2022-04-28-changing-spatial-resolution-of-a-raster-with-terra/
start_time=Sys.time()
deepwatersrs2=rast(nrows=DIMY, ncols=DIMX, extent=ext(deepwaterscrop2), crs=crs(deepwaterscrop2))
deepwatersrs2=resample(deepwaterscrop2, deepwatersrs2, method='bilinear')
end_time=Sys.time()
lapso_terra=end_time-start_time  # Time difference of 1.710264 secs
plot(deepwatersrs2)
abline(v=0)  # Greenwich meridian


# Improvement raster vs terra packages:

# action    raster      terra       unit    terra/raster
# --------- ----------- ----------- ------- ------------
# REPROJECT 4.666385    1.858117    min     39.8%
# PLOT      103.1214    10.92895    s       10.6%
# CROP      4.796154    2.65933     s       55.4%
# RESAMPLE  25.8529     1.710264    s       6.6%