# ------------------------------------------------------------------------------
# ---------- SET PARAMETERS  ---------------------------------------------------
# ------------------------------------------------------------------------------

# command line arguments
args = commandArgs(trailingOnly=TRUE)
for(i in seq_along(args)){
    if(args[i] == "-ngif") ngif = as.integer(args[i+1])
    if(args[i] == "-res") res = as.integer(args[i+1])
}
if(!exists(deparse(substitute(ngif)))) ngif = 10
if(!exists(deparse(substitute(res)))) res = 720

# ------------------------------------------------------------------------------
# ---------- LOAD THE RESULTS  -------------------------------------------------
# ------------------------------------------------------------------------------

# load the libraries
suppressMessages(suppressWarnings(library(INLA)))
library(ggplot2)
library(maps)
library(mapdata)

# load the data
load("data/data.Rdata")
mu.inv = readBin("../../data/mu", "double", n.s*n.t+n.b+1, endian = "swap")[-1]
ngif = min(ngif, n.t)

# ------------------------------------------------------------------------------
# ---------- CREATE PROJECTIONS  -----------------------------------------------
# ------------------------------------------------------------------------------

# create a US border
usa = map_data("usa")
usa = usa[1:6886,]
coordinates(usa) = ~ long + lat
usa@proj4string = CRS("+proj=longlat +datum=WGS84 +lon_0=100")
usa = spTransform(usa, CRS("+proj=moll +units=km"))
borderline = usa@coords
borderline[,1] = (borderline[,1]-min(bound$loc[,1])) / (diff(range(bound$loc[,1])))
borderline[,2] = (borderline[,2]-min(bound$loc[,2])) / (diff(range(bound$loc[,2])))

# create a projection
grid.ratio = diff(range(bound$loc[,2])) / diff(range(bound$loc[,1]))
grid.res = c(res, res)
grid.res[2] = round(grid.res[1] * grid.ratio)
grid.x = seq(min(bound$loc[,1]), max(bound$loc[,1]), length = grid.res[1])
grid.y = seq(min(bound$loc[,2]), max(bound$loc[,2]), length = grid.res[2])
grid.loc = expand.grid(grid.x, grid.y)
grid.loc = as.matrix(grid.loc)
grid.A = inla.spde.make.A(mesh = mesh.s, loc = grid.loc)

# set the plot range
zlim.mu = range(mu.inv[1:(n.s*n.t)])

# ------------------------------------------------------------------------------
# ---------- PLOT THE RESULTS  -------------------------------------------------
# ------------------------------------------------------------------------------

# choose time points
tt = 1:ngif

# plot several slices
for(i in seq_along(tt)){
    # set time
    t = tt[i]
    
    # project the mean
    mu.inv.grid = drop(grid.A %*% mu.inv[1:n.s + (t-1)*n.s])
    mu.inv.grid = matrix(mu.inv.grid, grid.res[1])
    
    # plot the spatial field
    png(paste0("img/",i,".png"), width=grid.res[1], height=grid.res[2])
    par(mar=c(0,0,0,0))
    image(mu.inv.grid, col=viridis::viridis(100), asp=grid.ratio, zlim=zlim.mu, xaxt="n", yaxt="n", axes=FALSE)
    lines(borderline, col="white", lwd=2)
    day = as.character.Date(as.Date(x = t, origin = "2010-12-31"), "%d %b %Y")
    text(x = 0.15, y = 0.05, day, col = "white", cex = 2)
    # mtext(day, side=1, line=0.5, at=0, cex=2)
    dev.off()
}

# generate a gif
system(paste("python3 plotgif.py ", ngif, "img", "field.gif"))

# clean up
for(i in seq_along(tt)){
    file.remove(paste0("img/",i,".png"))
}
