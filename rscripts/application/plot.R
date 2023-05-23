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
d.inv = readBin("../../data/d", "double", n.s*n.t+n.b+1, endian = "swap")[-1]
sd.inv = sqrt(d.inv)
save(list = c("mu.inv", "sd.inv"), file = "data/result.Rdata")

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
grid.res = c(256, 256)
grid.res[2] = round(grid.res[1] * grid.ratio)
grid.x = seq(min(bound$loc[,1]), max(bound$loc[,1]), length = grid.res[1])
grid.y = seq(min(bound$loc[,2]), max(bound$loc[,2]), length = grid.res[2])
grid.loc = expand.grid(grid.x, grid.y)
grid.loc = as.matrix(grid.loc)
grid.A = inla.spde.make.A(mesh = mesh.s, loc = grid.loc)

# set the plot range
zlim.mu = range(mu.inv[1:(n.s*n.t)])
zlim.sd = range(sd.inv[1:(n.s*n.t)])

# ------------------------------------------------------------------------------
# ---------- PLOT THE RESULTS  -------------------------------------------------
# ------------------------------------------------------------------------------

# choose time points
tt = c(1,2,4200,4201)

pdf("img/fig.app.1.pdf", width=14, height=5)
par(mfcol=c(2,length(tt)), mar=c(2,2,0,0))

# plot several slices
for(i in seq_along(tt)){
    # set time
    t = tt[i]
    
    # project the mean
    mu.inv.grid = drop(grid.A %*% mu.inv[1:n.s + (t-1)*n.s])
    mu.inv.grid = matrix(mu.inv.grid, grid.res[1])
    
    # project the sd
    sd.inv.grid = drop(grid.A %*% sd.inv[1:n.s + (t-1)*n.s])
    sd.inv.grid = matrix(sd.inv.grid, grid.res[1])
    
    # plot the spatial field
    image(mu.inv.grid, col=viridis::viridis(100), asp=grid.ratio, zlim=zlim.mu, xaxt="n", yaxt="n", axes=FALSE)
    lines(borderline, col="white", lwd=2)
    if(i==1) title(ylab="Mean", line=1)
    image(sd.inv.grid, col=viridis::inferno(100), asp=grid.ratio, zlim=zlim.sd, xaxt="n", yaxt="n", axes=FALSE)
    lines(borderline, col="white", lwd=2)
    if(i==1) title(ylab="SD", line=1)
    title(xlab=paste0("t = ", t), line = 1)
}
dev.off()
