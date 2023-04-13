source("plot.R")

bbound = bound
bbound$loc[,1] = (bbound$loc[,1]-min(range(lattice.long))) / diff(range(lattice.long))
bbound$loc[,2] = (bbound$loc[,2]-min(range(lattice.lat))) / diff(range(lattice.lat))
for(i in 0:364){
    mu.lattice = drop(A.lattice %*% mu[(1:n.s)+(n.s*i)])
    mu.lattice = matrix(mu.lattice, n.long)
    png(paste0("img/", i, ".png"), width = 1080, height = 720)
    image(mu.lattice, asp = lattice.ratio, col = viridis::viridis(200), zlim = range(mu))
    lines(bbound$loc, col="white", lwd=2)
    lines(bbound$loc[c(dim(bbound$loc)[1],1),], col="white", lwd=2)
    text(x = 0.95, y = 0.95, month.name[1+i%/%(365/12)], cex = 2)
    # Sys.sleep(0.05)
    dev.off()
}