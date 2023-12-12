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
    
    return(hillshade^(1/gamma))
}

NewBitmap = function(dimx, dimy, val=0) {
    # Crea bitmap de dimensiones dimx y dimy
    return(array(val,c(dimx,dimy)))
}

# Por Carlos Gil Bellosta
indices.drawline = function(x0, y0, x1, y1) {
    x0=round(x0)
    x1=round(x1)
    y0=round(y0)
    y1=round(y1)
    
    if (y0 == y1) return(cbind(x0:x1, y0)) # Recta de m=0 o un punto
    if (abs(x1 - x0) >= abs(y1 - y0)) { # Recta de 0 < |m| <= 1
        m = (y1 - y0) / (x1 - x0)
        cbind(x0:x1, round(y0 + m * ((x0:x1) - x0)))
    } else indices.drawline(y0, x0, y1, x1)[, 2:1]  # Recta de |m| > 1
    # Llamada traspuesta recursiva y traspuesta
}

DrawLine = function(img, x0, y0, x1, y1, inc=TRUE, val=1) {
    # Dibuja recta desde (x0,y0)-(x1,y1)
    # Por defecto método no destructivo y con valor=1
    indices=indices.drawline(x0, y0, x1, y1)
    if (inc) img[indices]=img[indices]+val
    else img[indices]=val
    
    return(img)
}

DrawRect = function(img, x0, y0, x1, y1, inc=TRUE, val=1, fill=FALSE) {
    # Dibuja rectángulo (x0,y0)-(x1,y1)
    # Por defecto método no destructivo, con valor=1 y sin relleno
    x0=round(x0)
    x1=round(x1)
    y0=round(y0)
    y1=round(y1)
    
    if (fill) {
        if (inc) img[x0:x1,y0:y1]=img[x0:x1,y0:y1]+val
        else img[x0:x1,y0:y1]=val
        
        return(img)
    } else {
        indices=which( ( (row(img)==x0         | row(img)==x1        ) &
                         (col(img)>=min(y0,y1) & col(img)<=max(y0,y1)) ) |
                       ( (col(img)==y0         | col(img)==y1        ) &
                         (row(img)>=min(x0,x1) & row(img)<=max(x0,x1)) ) )
        if (inc) img[indices]=img[indices]+val
        else img[indices]=val
        
        return(img)
    }
}

DibujarNumero = function(img, x0, y0, inc=FALSE, val=1, fill=FALSE,
                         num, width, height) {
    # Dibuja cifra 0-9 en (x0,y0)
    # Por defecto método no destructivo y con valor=1
    
    if (num=='0') { 
        img=DrawRect(img, x0, y0, x0+width, y0-height, inc, val, fill)
    } else if (num=='1') {
        img=DrawLine(img, x0+width/2, y0, x0+width/2, y0-height, inc, val)
    } else if (num=='2') {
        img=DrawLine(img, x0, y0, x0+width, y0, inc, val)
        img=DrawLine(img, x0+width, y0, x0+width, y0-height/2, inc, val)
        img=DrawLine(img, x0+width, y0-height/2, x0, y0-height/2, inc, val)
        img=DrawLine(img, x0, y0-height/2, x0, y0-height, inc, val)
        img=DrawLine(img, x0, y0-height, x0+width, y0-height, inc, val)
    } else if (num=='3') {
        img=DrawLine(img, x0, y0, x0+width, y0, inc, val)
        img=DrawLine(img, x0, y0-height/2, x0+width, y0-height/2, inc, val)
        img=DrawLine(img, x0, y0-height, x0+width, y0-height, inc, val)
        img=DrawLine(img, x0+width, y0, x0+width, y0-height, inc, val)
    } else if (num=='4') {
        img=DrawLine(img, x0, y0, x0, y0-height/2, inc, val)
        img=DrawLine(img, x0, y0-height/2, x0+width, y0-height/2, inc, val)
        img=DrawLine(img, x0+width, y0, x0+width, y0-height, inc, val)
    } else if (num=='5') {
        img=DrawLine(img, x0+width, y0, x0, y0, inc, val)
        img=DrawLine(img, x0, y0, x0, y0-height/2, inc, val)
        img=DrawLine(img, x0, y0-height/2, x0+width, y0-height/2, inc, val)
        img=DrawLine(img, x0+width, y0-height/2, x0+width, y0-height, inc, val)
        img=DrawLine(img, x0+width, y0-height, x0, y0-height, inc, val)
    } else if (num=='6') {
        img=DrawRect(img, x0, y0-height/2, x0+width, y0-height, inc, val, fill)
        img=DrawLine(img, x0, y0, x0+width, y0, inc, val)
        img=DrawLine(img, x0, y0, x0, y0-height/2, inc, val)
    } else if (num=='7') {
        img=DrawLine(img, x0, y0, x0+width, y0, inc, val)
        img=DrawLine(img, x0+width, y0, x0+width, y0-height, inc, val)
    } else if (num=='8') {
        img=DrawRect(img, x0, y0, x0+width, y0-height/2, inc, val, fill)
        img=DrawRect(img, x0, y0-height/2, x0+width, y0-height, inc, val, fill)
    } else if (num=='9') {
        img=DrawRect(img, x0, y0, x0+width, y0-height/2, inc, val, fill)
        img=DrawLine(img, x0+width, y0-height/2, x0+width, y0-height, inc, val)
        img=DrawLine(img, x0, y0-height, x0+width, y0-height, inc, val)
    } else if (num=='-') {
        img=DrawLine(img, x0, y0-height/2, x0+width, y0-height/2, inc, val)
    } else if (num=='m') {
        img=DrawLine(img, x0, y0-height/2, x0+width, y0-height/2, inc, val)
        img=DrawLine(img, x0, y0-height/2, x0, y0-height, inc, val)
        img=DrawLine(img, x0+width/2, y0-height/2, x0+width/2, y0-height, inc, val)
        img=DrawLine(img, x0+width, y0-height/2, x0+width, y0-height, inc, val)
    } else {
        return(img)  # Cifra inválida
    }

    return(img)
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
# outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)]+
#     outline[1:(DIMY-2), 2:(DIMX-1)]+outline[3:(DIMY-0), 2:(DIMX-1)]+
#     outline[2:(DIMY-1), 1:(DIMX-2)]+outline[2:(DIMY-1), 3:(DIMX-0)]

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
COLUP=1.8  # 2
COLDOWN=0.56  # 0.5
WATERMIN=0.3  # difusse under water colours a bit
WATERMAX=0.7
BORDERCOLOUR=0.2

mappre=array(0, c(DIMY, DIMX, 3))
mappre[,,1]=hill
mappre[,,2]=hill
mappre[,,3]=hill
border=which(outline==1)  # to draw map borders

# Initial sequence (water spreads)
NFRAMES=2190  # 73.00s audio track at 30fps
RATIO=2  # to slow down initial phases
for (frame in 0:(NFRAMES-1)) {
    WATERLEVEL=MAXHEIGHT * (frame/(NFRAMES-1))^RATIO

    land=which(DEM>=WATERLEVEL)
    water=which(DEM<WATERLEVEL)
    map=mappre
    map[,,1][land]=map[,,1][land]^(1/COLUP)  # R land
    map[,,3][land]=map[,,3][land]^(1/COLDOWN)  # B land
    
    map[,,1][water]=map[,,1][water]^(1/COLDOWN)  # R water
    map[,,1][water]=(WATERMAX-WATERMIN)*map[,,1][water]+WATERMIN
    map[,,2][water]=(WATERMAX-WATERMIN)*map[,,2][water]+WATERMIN
    map[,,3][water]=(WATERMAX-WATERMIN)*map[,,3][water]+WATERMIN
    
    map[,,1][border]=BORDERCOLOUR
    map[,,2][border]=BORDERCOLOUR
    map[,,3][border]=BORDERCOLOUR
    
    # Write water level label
    TXT=paste0(as.character(round(WATERLEVEL)),'m')
    LONG=nchar(TXT)
    label=NewBitmap(251, 63)
    for (i in 1:LONG) {
        num=substring(TXT, i, i)
        label=DibujarNumero(label, 2+i*44-44, 62, num=num, width=28, height=60)
    }
    label=t(label[,ncol(label):1])
    # Stroke=3
    label[1:61,2:250]=label[1:61,2:250]+label[2:62,2:250]
    label[3:63,2:250]=label[3:63,2:250]+label[2:62,2:250]
    label[2:62,1:249]=label[2:62,1:249]+label[2:62,2:250]
    label[2:62,3:251]=label[2:62,3:251]+label[2:62,2:250]
    meters=which(label!=0)
    map[950:(950+63-1),(946-(LONG-1)/2*44):(946-(LONG-1)/2*44+251-1),1][meters]=0
    map[950:(950+63-1),(946-(LONG-1)/2*44):(946-(LONG-1)/2*44+251-1),2][meters]=0
    map[950:(950+63-1),(946-(LONG-1)/2*44):(946-(LONG-1)/2*44+251-1),3][meters]=0
    
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
    land=which(DEM>=WATERLEVEL)
    water=which(DEM<WATERLEVEL)
    map=mappre
    map[,,1][land]=map[,,1][land]^(1/COLUP)  # R land
    map[,,3][land]=map[,,3][land]^(1/COLDOWN)  # B land
    
    map[,,1][water]=map[,,1][water]^(1/COLDOWN)  # R water
    map[,,1][water]=(WATERMAX-WATERMIN)*map[,,1][water]+WATERMIN
    map[,,2][water]=(WATERMAX-WATERMIN)*map[,,2][water]+WATERMIN
    map[,,3][water]=(WATERMAX-WATERMIN)*map[,,3][water]+WATERMIN

    map[,,1][border]=BORDERCOLOUR
    map[,,2][border]=BORDERCOLOUR
    map[,,3][border]=BORDERCOLOUR
    
    # Write water level label
    TXT=paste0(as.character(round(WATERLEVEL)),'m')
    LONG=nchar(TXT)
    label=NewBitmap(251, 63)
    for (i in 1:LONG) {
        num=substring(TXT, i, i)
        label=DibujarNumero(label, 2+i*44-44, 62, num=num, width=28, height=60)
    }
    label=t(label[,ncol(label):1])
    # Stroke=3
    label[1:61,2:250]=label[1:61,2:250]+label[2:62,2:250]
    label[3:63,2:250]=label[3:63,2:250]+label[2:62,2:250]
    label[2:62,1:249]=label[2:62,1:249]+label[2:62,2:250]
    label[2:62,3:251]=label[2:62,3:251]+label[2:62,2:250]
    meters=which(label!=0)
    map[950:(950+63-1),(946-(LONG-1)/2*44):(946-(LONG-1)/2*44+251-1),1][meters]=0
    map[950:(950+63-1),(946-(LONG-1)/2*44):(946-(LONG-1)/2*44+251-1),2][meters]=0
    map[950:(950+63-1),(946-(LONG-1)/2*44):(946-(LONG-1)/2*44+251-1),3][meters]=0

    writePNG(map, paste0("img", ifelse(frame+Offset<10, "000",
                ifelse(frame+Offset<100, "00",
                ifelse(frame+Offset<1000, "0", ""))), frame+Offset, ".png"))
    
    print(paste0(frame+1, "/", NFRAMES, ": Water level ", WATERLEVEL))
}



# Building MP4:
# ffmpeg -framerate 30 -i img%4d.png -i oblivionwakingup.wav
#        -c:v libx264 -crf 18 -pix_fmt yuv420p deepwaters.mp4




