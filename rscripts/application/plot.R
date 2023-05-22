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
load("data/theta.Rdata")
load("data/result.Rdata")
if(TRUE){           # rewrite?
    mu.inv = readBin("../../data/mu", "double", n.s*n.t+n.b+1, endian = "swap")[-1]
    d.inv = readBin("../../data/d", "double", n.s*n.t+n.b+1, endian = "swap")[-1]
    sd.inv = sqrt(d.inv)
    save(list = c("mu.inla", "sd.inla", "mu.inv", "sd.inv"), file = "data/result.Rdata")
} else {
    mu.inv = mu.inla
    sd.inv = sd.inla
}

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
zlim.mu = range(mu.inla[1:(n.s*n.t)], mu.inv[1:(n.s*n.t)])
zlim.sd = range(sd.inla[1:(n.s*n.t)], sd.inv[1:(n.s*n.t)])

# ------------------------------------------------------------------------------
# ---------- PLOT THE RESULTS  -------------------------------------------------
# ------------------------------------------------------------------------------

# choose time points
tt = round(seq(1,n.t,length=10))

# plot several slices
for(i in seq_along(tt)){
    # set time
    t = tt[i]
    
    # project the mean
    mu.inla.grid = drop(grid.A %*% mu.inla[1:n.s + (t-1)*n.s])
    mu.inla.grid = matrix(mu.inla.grid, grid.res[1])
    mu.inv.grid = drop(grid.A %*% mu.inv[1:n.s + (t-1)*n.s])
    mu.inv.grid = matrix(mu.inv.grid, grid.res[1])
    
    # project the sd
    sd.inla.grid = drop(grid.A %*% sd.inla[1:n.s + (t-1)*n.s])
    sd.inla.grid = matrix(sd.inla.grid, grid.res[1])
    sd.inv.grid = drop(grid.A %*% sd.inv[1:n.s + (t-1)*n.s])
    sd.inv.grid = matrix(sd.inv.grid, grid.res[1])
    
    # plot the spatial field
    # png(paste0("img/fig.app.field.",i,".png"), width=1024, height=1024)
    pdf(paste0("img/fig.app.field.",i,".pdf"), width=5, height=5)
    par(mfrow=c(2,2), mar=c(3,3,1,1))
    image(mu.inla.grid, col=viridis::viridis(100), zlim=zlim.mu, xaxt="n", yaxt="n")
    lines(borderline, col="white", lwd=2)
    title(ylab="Mean", line=1)
    image(mu.inv.grid, col=viridis::viridis(100), zlim=zlim.mu, xaxt="n", yaxt="n")
    lines(borderline, col="white", lwd=2)
    image(sd.inla.grid, col=viridis::inferno(100), zlim=zlim.sd, xaxt="n", yaxt="n")
    lines(borderline, col="white", lwd=2)
    title(xlab="R-INLA", ylab="SD", line=1)
    image(sd.inv.grid, col=viridis::inferno(100), zlim=zlim.sd, xaxt="n", yaxt="n")
    lines(borderline, col="white", lwd=2)
    title(xlab="Ovelapping RBMC", line=1)
    dev.off()
}

# plot the temporal field
pdf("img/fig.app.2.pdf", width=5, height=5)
par(mfrow=c(1,1), mar=c(2,2,1,1))
t.idx = (1:n.t-1)*n.s+111 # +111 = spatial point 111
error = abs(sd.inv/sd.inla-1)[t.idx]
plot(error, t="l", xlab=NA, ylab=NA, xaxt="n", yaxt="n", ylim=c(0,max(error)))
axis(1, at=seq(0,n.t,length=5), labels=seq(0,n.t,length=5)*c(1,1,NA,1,1))
axis(2, at=seq(0,max(error),length=2), labels=c(0, formatC(max(error), digits=2, format="e")))
title(xlab="Time", ylab="Relative error", line=1)
dev.off()