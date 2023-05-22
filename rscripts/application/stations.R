### consider US stations information 

library(sp)

if(!any(ls()=="dpath"))
    dpath <- "data"

if(file.exists(
    file.path(dpath, "stations.rds"))) {

    stations <- readRDS(
        file.path(dpath, "stations.rds"))

} else {

    if(file.exists(
        file.path(dpath, "allstations.rds"))) {

        allstations <- readRDS(
            file.path(dpath, "allstations.rds"))

    } else {
### stations info from
### ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt
        
### check if the file is available locally
        stfl <- 'ghcnd-stations.txt'
        if(!file.exists(file.path(dpath, stfl)))
            download.file(paste0('ftp://ftp.ncdc.noaa.gov/',
                                 'pub/data/ghcn/daily/', stfl),
                          file.path(dpath, stfl))
        
### width of the colums in the file:          
        (ws <- diff(c(0,11,20,30,37,71,75,79,85)))
        
### read station information: longitute, latitude & altitude information
        allstations <- read.fwf(
            file.path(dpath, stfl), ws[1:4])
        colnames(allstations) <- c('station', 'latitude', 'longitude', 'elevation')
    
        dim(allstations)

### deal with sp and projection
        coordinates(allstations) <- ~ longitude + latitude
        allstations@proj4string <- CRS('+proj=longlat +datum=WGS84 +lon_0=100')

        head(allstations)

        saveRDS(allstations,
                file.path(dpath, "allstations.rds"))

    }
    
### index of stations with 'US' code
    ii0us <- which(substr(allstations$station, 1, 2)=='US')
    table(substr(allstations$station, 1, 3)[ii0us])
    
    ## define a box around US main territory
    box0ll <- SpatialPolygons(list(Polygons(list(Polygon(
        cbind(c(-130, -60, -60, -130, -130),
              c(50, 50, 23, 23, 50)))), '0')),
        proj4string=allstations@proj4string)
    
    ii1us <- which(!is.na(over(allstations[ii0us,], box0ll)))
    str(ii1us)
    
    if(FALSE) 
        plot(allstations[ii0us[ii1us],])
    
### projection with units in km
    stations <- spTransform(
        allstations[ii0us[ii1us], ],
        CRS('+proj=moll +units=km'))
    
    if(FALSE) {
        
        par(mfrow=c(1,1), mar=c(0,0,0,0))
        plot(allstations, cex=0.3, pch=19)
        
    }
    
    saveRDS(stations,
            file.path(dpath, "stations.rds"))
}

