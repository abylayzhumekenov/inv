# fit the model
m.s = 1000
m.t = 365
n.cores1 = 2
n.cores2 = 6
source("fitinla.R")

# read results
mu = result$summary.random$field$mode
# mu = readBin("../../data/mu", "double", 1+n.st+n.b, endian="swap")[-1]
# d = readBin("../../data/d", "double", 1+n.st+n.b, endian="swap")[-1]

# create a grid
grid.ratio = diff(range(mesh.s$loc[,2])) / diff(range(mesh.s$loc[,1]))
grid.nx = 1000
grid.ny = round(grid.nx * grid.ratio)
grid.x = seq(min(mesh.s$loc[,1]), max(mesh.s$loc[,1]), length = grid.nx)
grid.y = seq(min(mesh.s$loc[,2]), max(mesh.s$loc[,2]), length = grid.ny)
grid.loc = expand.grid(grid.x, grid.y)
# colnames(grid.loc) = c("xcoord", "ycoord")
grid.loc = as.matrix(grid.loc)
grid.A = inla.spde.make.A(mesh = mesh.s, loc = grid.loc)

# project solution
grid.mu = drop(grid.A %*% mu[(1:n.s)+(n.s*364)])
grid.mu = matrix(grid.mu, grid.nx)
image(grid.mu, asp = grid.ratio, col = viridis::viridis(100))

# d.lattice = drop(A.lattice %*% d[(1:n.s)+(n.s*1000)])
# d.lattice = matrix(d.lattice, n.long)
# image(d.lattice, asp = lattice.ratio, col = rev(viridis::inferno(100)))

# generate gif
gif.plot = TRUE
if(gif.plot){
    gif.width = 720
    gif.height = 480
    source("plotgif.R")
    system(paste0("python3 plotgif.py ", n.t))
}
