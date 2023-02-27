
### years of data needed
fnyears <- function(n) {
    Dates <- as.Date('2022-12-31')-(n-1):0
    unique(substr(as.character(Dates), 1, 4))
}

### download data (done once)
fdownload <- function(year, dpath) {
    if(!file.exists(file.path(
            dpath, paste0(year, '.csv.gz'))))
        download.file(
            paste0("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/", 
                   "daily/by_year/", year, '.csv.gz'),
            file.path(dpath, paste0(year, '.csv.gz')))
}

### get the long format data for each year and reshape 
wdatf <- function(year, dpath, vnam, verbose=TRUE) {
    if(length(vnam)>0) {
        vnam <- vnam[1]
        warning("length(vnam)>1 and taken vnam[1]!")
    } 
    wdfl <- file.path(dpath, paste0(vnam, 'wd', year, 'US.RData'))
    if(!file.exists(wdfl)) {
        fdownload(year, dpath)
        t0 <- Sys.time()
        if(verbose) 
            cat('reading data from year', year, '\n')
### read the data
        dat <- as.data.frame(fread(
            file.path(dpath, paste0(year, '.csv.gz')),
            nThread=1L))
        t1 <- Sys.time()
        if(verbose) {
            cat('year', year, 'readed', nrow(dat), 'lines ')
            print(t1-t0)
        }
### index of US stations with V6='' (quality control)
        i0sel <- (dat$V6=='') &
            (dat$V1%in%stations$station)
### and 'vnam'
        iidsel <- which(i0sel & (dat$V3 %in% vnam))
        t2 <- Sys.time()
        if(verbose) {
            cat('year', year, 'selected', length(iidsel), 'obs. ')
            print(t2-t1)
        }
        dat <- reshape(
            dat[iidsel, c(1,2,4)], dir='wide',
            timevar='V2', idvar='V1')
        colnames(dat) <- gsub(
            'V4.', '', colnames(dat), fixed=TRUE)
        t3 <- Sys.time()
        if(verbose) {
            cat('year', year, 'reshape', dim(dat[, -1]), 'obs. ')
            print(t3-t2)
        }
        if(is.character(dat[,1])) {
            rownames(dat) <- dat[,1]
            dat <- as.matrix(dat[,-1])
        }
        save(list='dat', file=wdfl)
        t4 <- Sys.time()
        if(verbose) {
            cat(wdfl, 'saved ')
            print(t4-t3)
        }
    }
    print(system(paste('ls -lh', wdfl)))
}
