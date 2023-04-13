dir.create("img")

boundline = bound$loc
boundline[,1] = (boundline[,1]-min(grid.x)) / diff(range(grid.x))
boundline[,2] = (boundline[,2]-min(grid.y)) / diff(range(grid.y))
boundline = rbind(boundline, boundline[1,])
for(i in 1:n.t){
    grid.mu = drop(grid.A %*% mu[(1:n.s)+(n.s*(i-1))])
    grid.mu = matrix(grid.mu, grid.nx)
    
    png(paste0("img/", i, ".png"), width = gif.width, height = gif.height)
    image(grid.mu, asp = grid.ratio, col = viridis::viridis(200), zlim = range(mu))
    lines(boundline, col = "white", lwd = 2)
    day = as.character.Date(as.Date(x = i, origin = "2010-12-31"), "%d %b %Y")
    text(x = 0.85, y = 0.05, day, col = "white", cex = 2)
    dev.off()
}