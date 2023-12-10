# Deep Waters. Practice with GeoTIFF format and raster package
# www.overfitting.net
# https://www.overfitting.net/

    
library(raster)  # read GeoTIFF, reprojection, crop and resample
library(tiff)  # save 16-bit TIFF's
library(png)  # save 8-bit PNG's


hillshade=function(map, dx=25, dlight=c(0, 2, 3), gamma=1) {
    # dx: map resolution (m)
    # Lighting direction (3D vector defined from observer to light source):
    #   dlight=c(0, 2, 3)  # sunrise
    #   dlight=c(0, 0, 1)  # midday
    #   dlight=c(0,-2, 3)  # sunset
    
    DIMY=nrow(map)    
    DIMX=ncol(map)
    
    # Save map
    mapnorm=map-min(map)
    writeTIFF((mapnorm/max(mapnorm))^(1/gamma), "map.tif",
              bits.per.sample=16, compression="LZW")
    rm(mapnorm)
    
    dlightM=sum(dlight^2)^0.5
    
    # Vectorial product to calculate n (orthogonal vector)
    nx = 2*dx*(map[1:(DIMY-2), 2:(DIMX-1)] - map[3:DIMY,     2:(DIMX-1)])
    ny = 2*dx*(map[2:(DIMY-1), 1:(DIMX-2)] - map[2:(DIMY-1), 3:DIMX])
    nz = 4*dx^2
    nM = (nx^2 + ny^2 + nz^2)^0.5
    
    # Dot product to calculate cos(theta)
    dn = dlight[1]*nx + dlight[2]*ny + dlight[3]*nz  # (DIMY-2)x(DIMX-2) matrix
    
    # Reflectance (=cos(theta))
    hillshadepre=dn/(dlightM*nM)
    hillshadepre[hillshadepre<0]=0  # clip negative values
    
    # Add 'lost' borders
    hillshade=matrix(0, nrow=DIMY, ncol=DIMX)
    hillshade[2:(DIMY-1), 2:(DIMX-1)]=hillshadepre
    rm(hillshadepre)
    hillshade[c(1,DIMY),]=hillshade[c(2,DIMY-1),]
    hillshade[,c(1,DIMX)]=hillshade[,c(2,DIMX-1)]
    
    return(hillshade)
}


###########################################################

# 1. PROCESS GEOTIFF DATA TO GET THE DEM AS A FULL HD RESOLUTION MATRIX

# https://download.gebco.net/
# The GEBCO_2023 Grid is a global terrain model for ocean and land,
# providing elevation data, in meters, on a 15 arc-second interval grid
# of 43200 rows x 86400 columns, giving 3,732,480,000 data points.
# The data values are pixel-centre registered i.e. they refer to elevations,
# in meters, at the centre of grid cells.
deepwaters=raster("gebco_2023_n69.5654_s25.1807_w-25.752_e31.2891.tif")  # read GeoTIFF file
deepwaters
plot(deepwaters)
abline(v=0)  # Greenwich meridian


# REPROJECT raster from Longitude Latitude (+proj=longlat)/WGS84
# to Lambert Conic Conformal (+proj=lcc)/WGS84
# https://pygis.io/docs/d_understand_crs_codes.html
# https://stackoverflow.com/questions/36868506/how-to-change-a-lambert-conic-conformal-raster-projection-to-latlon-degree-r
CRS="+proj=lcc +units=km"  # by default crs="+proj=lcc +ellps=WGS84 +lat_1=33 +lat_2=45 +lon_0=0"
deepwatersrp=projectRaster(deepwaters, crs=CRS)
deepwatersrp
plot(deepwatersrp)
abline(v=0)  # Greenwich meridian


# CROP raster to area of interest (save Canary Islands and The Netherlands)
cropdef=as(extent(-1850, 2200, 3450, 6050), 'SpatialPolygons')  # kms
crs(cropdef)=crs(deepwatersrp)
deepwaterscrop=crop(deepwatersrp, cropdef)
deepwaterscrop
plot(deepwaterscrop)
abline(v=0)  # Greenwich meridian


# RESAMPLE raster to Full HD
DIMY=1080
DIMX=1920
deepwatersrs=raster(nrow=DIMY, ncol=DIMX,
                    ext=extent(deepwaterscrop), crs=crs(deepwaterscrop))
deepwatersrs=resample(deepwaterscrop, deepwatersrs,
                      method='bilinear')  # bilinear is faster than ngb!
plot(deepwatersrs)
abline(v=0)  # Greenwich meridian


# Convert to matrix and save as TIFF
DEM=as.matrix(deepwatersrs)
hist(DEM, breaks=1000)
abline(v=0, col='red')
# DEM[is.na(DEM)]=0  # no NA's
DEM=DEM-min(DEM)
writeTIFF(DEM/max(DEM), "deepwaters.tif",
          compression='LZW', bits.per.sample=16)

DEM=as.matrix(deepwatersrs)
rm(deepwaters)
rm(deepwatersrp)
rm(deepwaterscrop)
rm(deepwatersrs)


###########################################################

# 2. PROCESS MATRIX TO OBTAIN MAP CONTOURS AND HILLSHADE

# Calculate solid map contour
solid=DEM
solid[solid>=0]=1  # set >=0 areas to 1 (land)
solid[solid<0]=0  # set <0 areas to 0 (water)
writePNG(solid, "mapsolid.png")

# Calculate outline map from solid map
outline=solid*0
# 1 pixel thickness outline
outline[2:(DIMY-1), 2:(DIMX-1)]=
    abs(solid[1:(DIMY-2), 2:(DIMX-1)] -
            solid[2:(DIMY-1), 2:(DIMX-1)]) +
    abs(solid[2:(DIMY-1), 1:(DIMX-2)] -
            solid[2:(DIMY-1), 2:(DIMX-1)])
# increase to 2 pixel thickness outline
outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)]+
    outline[1:(DIMY-2), 2:(DIMX-1)]+outline[2:(DIMY-1), 3:(DIMX-0)]
# increase to 3 pixel thickness outline
outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)]+
    outline[1:(DIMY-2), 2:(DIMX-1)]+outline[3:(DIMY-0), 2:(DIMX-1)]+
    outline[2:(DIMY-1), 1:(DIMX-2)]+outline[2:(DIMY-1), 3:(DIMX-0)]

outline[outline!=0]=1
writePNG(outline, "mapoutline.png")


# Calculate grayscale hillshade
MIX=0.85  # two light sources are mixed to fill dark areas a bit
hill=hillshade(DEM, dx=800, dlight=c(0, 2, 3))
hillfill=hillshade(DEM, dx=800, dlight=c(0, 3, 2))
hill=hill*MIX+hillfill*(1-MIX)
gamma=1/2.0
hill=(hill/max(hill))^(1/gamma)  # darken hillshade a bit

# Save hillshade
writeTIFF(hill, "hillshade.tif",
          bits.per.sample=16, compression="LZW")

# Display hillshade
image(t(hill[nrow(hill):1,]), useRaster=TRUE,
      col=c(gray.colors(256, start=0, end=1, gamma=0.5)),
      asp=nrow(hill)/ncol(hill), axes=FALSE)


###########################################################

# 3. GENERATE ANIMATION

MAXHEIGHT=max(DEM)
MINHEIGHT=min(DEM)
COLUP=2
COLDOWN=0.5
WATERMIN=0.2  # to difusse under water colours
WATERMAX=0.8

mappre=array(0, c(DIMY, DIMX, 3))
mappre[,,1]=hill
mappre[,,2]=hill
mappre[,,3]=hill
indicesborder=which(outline==1)  # draw map borders
mappre[,,1][indicesborder]=0
mappre[,,2][indicesborder]=0
mappre[,,3][indicesborder]=0

# Initial sequence (water spreads)
NFRAMES=2190  # 73.00s audio track at 30fps
RATIO=2  # to slow down initial phases
for (frame in 0:(NFRAMES-1)) {
    WATERLEVEL=MAXHEIGHT * (frame/(NFRAMES-1))^RATIO
    indicesland=which(DEM>=WATERLEVEL)
    indiceswater=which(DEM<WATERLEVEL)
    map=mappre
    map[,,1][indicesland]=map[,,1][indicesland]^(1/COLUP)  # R land
    map[,,3][indicesland]=map[,,3][indicesland]^(1/COLDOWN)  # B land
    
    map[,,1][indiceswater]=map[,,1][indiceswater]^(1/COLDOWN)  # R water
    map[,,1][indiceswater]=(WATERMAX-WATERMIN)*map[,,1][indiceswater]+WATERMIN
    map[,,2][indiceswater]=(WATERMAX-WATERMIN)*map[,,2][indiceswater]+WATERMIN
    map[,,3][indiceswater]=(WATERMAX-WATERMIN)*map[,,3][indiceswater]+WATERMIN
    
    writePNG(map, paste0("img", ifelse(frame<10, "000",
                ifelse(frame<100, "00",
                ifelse(frame<1000, "0", ""))), frame, ".png"))
    
    print(paste0(frame+1, "/", NFRAMES, ": Water level ", WATERLEVEL))
    
}

# Final sequence (land spreads)
Offset=2190  # previously generated frames
NFRAMES=544  # 18.13s audio track at 30fps
RATIO=1/2  # to speed up initial phases
for (frame in 0:(NFRAMES-1)) {
    WATERLEVEL=(MINHEIGHT-MAXHEIGHT) * (frame/(NFRAMES-1))^RATIO + MAXHEIGHT
    indicesland=which(DEM>=WATERLEVEL)
    indiceswater=which(DEM<WATERLEVEL)
    map=mappre
    map[,,1][indicesland]=map[,,1][indicesland]^(1/COLUP)  # R land
    map[,,3][indicesland]=map[,,3][indicesland]^(1/COLDOWN)  # B land
    
    map[,,1][indiceswater]=map[,,1][indiceswater]^(1/COLDOWN)  # R water
    map[,,1][indiceswater]=(WATERMAX-WATERMIN)*map[,,1][indiceswater]+WATERMIN
    map[,,2][indiceswater]=(WATERMAX-WATERMIN)*map[,,2][indiceswater]+WATERMIN
    map[,,3][indiceswater]=(WATERMAX-WATERMIN)*map[,,3][indiceswater]+WATERMIN
    
    writePNG(map, paste0("img", ifelse(frame+Offset<10, "000",
                ifelse(frame+Offset<100, "00",
                ifelse(frame+Offset<1000, "0", ""))), frame+Offset, ".png"))
    
    print(paste0(frame+1, "/", NFRAMES, ": Water level ", WATERLEVEL))
    
}


