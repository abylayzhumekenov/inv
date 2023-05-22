library(data.table)
library(parallel)
options(timeout = 1000)
ncores = 1

# command line arguments
args = commandArgs(trailingOnly=TRUE)
for(i in seq_along(args)){
	if(args[i] == "-ncores") ncores = as.integer(args[i+1])
}

source("functions.R")

if(!any(ls()=="dpath"))
    dpath <- "data"

if(!any(ls()=='ndays'))
    ndays <- 365*12+3

years <- fnyears(ndays)
names(years) <- years
years

source("stations.R")

### select the data
system.time(mclapply(
    sort(years), wdatf, dpath = dpath, vnam = "TMAX", verbose=TRUE,
    mc.cores=min(ncores, length(years), detectCores())))

### read each year data, select stations with at least 200 obs days
lwdat <- mclapply(sort(years), function(y) {
    wdfl <- file.path(dpath, paste0('TMAXwd', y, 'US.RData'))
    load(wdfl)    
    return(dat[which(rowSums(!is.na(dat))>200),])
}, mc.cores=min(ncores, length(years), detectCores()))

sapply(lwdat, dim)

### recursivelly merging (intersect)
wdat <- data.frame(station=rownames(lwdat[[1]]), lwdat[[1]])
if(length(lwdat)>1) {
    for(k in 2:length(lwdat)) {
        cat(sum(rownames(lwdat[[k]])%in%wdat[,1]),
            sum(wdat[,1]%in%rownames(lwdat[[k]])),
            length(intersect(wdat[,1], rownames(lwdat[[k]]))),
            length(union(wdat[,1], rownames(lwdat[[k]]))), '\n')        
        wdat <- merge(wdat, data.frame(station=rownames(lwdat[[k]]), lwdat[[k]]))
    }
}

rm(lwdat)
gc(reset=TRUE)

dim(wdat)
wdat[1:5, 1:5]

##plot(colSums(is.na(wdat[,-1])), type='l')
##plot(colSums(!is.na(wdat[,-1])), type='l')

### id to link data matrix to stations
ii.d.s <- pmatch(
    wdat[,1], 
    stations$station)
summary(ii.d.s)

stations.d <- stations[ii.d.s,]

if(FALSE)
    plot(stations.d)

cat("Matrix of data at", nrow(wdat), "stations on", ncol(wdat), "days\n")
