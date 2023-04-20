plot.gif = function(mu, bound, mesh.s, mesh.t, col, width = 720, height = 480, dirname = "img", filename = "field.gif"){
    n.s = mesh.s$n
    n.t = mesh.t$n

    # create a grid
    grid.ratio = diff(range(mesh.s$loc[,2])) / diff(range(mesh.s$loc[,1]))
    grid.nx = 1000
    grid.ny = round(grid.nx * grid.ratio)
    grid.x = seq(min(mesh.s$loc[,1]), max(mesh.s$loc[,1]), length = grid.nx)
    grid.y = seq(min(mesh.s$loc[,2]), max(mesh.s$loc[,2]), length = grid.ny)
    grid.loc = expand.grid(grid.x, grid.y)
    grid.loc = as.matrix(grid.loc)
    grid.A = inla.spde.make.A(mesh = mesh.s, loc = grid.loc)
    
    boundline = bound$loc
    boundline[,1] = (boundline[,1]-min(grid.x)) / diff(range(grid.x))
    boundline[,2] = (boundline[,2]-min(grid.y)) / diff(range(grid.y))
    boundline = rbind(boundline, boundline[1,])
    
    dir.create(dirname)
    for(i in 1:n.t){
        grid.mu = drop(grid.A %*% mu[(1:n.s)+(n.s*(i-1))])
        grid.mu = matrix(grid.mu, grid.nx)
        
        png(paste0(dirname, "/", i, ".png"), width = width, height = height)
        image(grid.mu, asp = grid.ratio, col = col, zlim = range(mu))
        lines(boundline, col = "white", lwd = 2)
        day = as.character.Date(as.Date(x = i, origin = "2010-12-31"), "%d %b %Y")
        text(x = 0.85, y = 0.05, day, col = "white", cex = 2)
        dev.off()
    }
    system(paste("python3 plotgif.py ", n.t, dirname, filename))
}
