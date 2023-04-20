# fit the model
m.s = 1000
m.t = 365
n.cores1 = 2
n.cores2 = 6
source("generate.R")
source("fitinla.R")
source("plotgif.R")

# read results
n.st = n.s * n.t
n.b = dim(data)[2] - 3
mu.inla = result$summary.random$field$mode
d.inla = result$summary.random$field$sd^2
sd.inla = sqrt(d.inla)
mu.inv = readBin("../../data/mu", "double", 1+n.st+n.b, endian="swap")[-1][1:n.st]
d.inv = readBin("../../data/d", "double", 1+n.st+n.b, endian="swap")[-1][1:n.st]
sd.inv = sqrt(d.inv)

# create a grid
grid.ratio = diff(range(mesh.s$loc[,2])) / diff(range(mesh.s$loc[,1]))
grid.nx = 1000
grid.ny = round(grid.nx * grid.ratio)
grid.x = seq(min(mesh.s$loc[,1]), max(mesh.s$loc[,1]), length = grid.nx)
grid.y = seq(min(mesh.s$loc[,2]), max(mesh.s$loc[,2]), length = grid.ny)
grid.loc = expand.grid(grid.x, grid.y)
grid.loc = as.matrix(grid.loc)
grid.A = inla.spde.make.A(mesh = mesh.s, loc = grid.loc)

# project solutions
t = 0

grid.mu.inla = drop(grid.A %*% mu.inla[(1:n.s)+(n.s*t)])
grid.mu.inla = matrix(grid.mu.inla, grid.nx)
image(grid.mu.inla, asp = grid.ratio, col = viridis::viridis(100), zlim = range(mu.inla))

grid.mu.inv = drop(grid.A %*% mu.inv[(1:n.s)+(n.s*t)])
grid.mu.inv = matrix(grid.mu.inv, grid.nx)
image(grid.mu.inv, asp = grid.ratio, col = viridis::viridis(100), zlim = range(mu.inv))

grid.sd.inla = drop(grid.A %*% sd.inv[(1:n.s)+(n.s*t)])
grid.sd.inla = matrix(grid.sd.inla, grid.nx)
image(grid.sd.inla, asp = grid.ratio, col = viridis::inferno(100), zlim = range(sd.inla))

grid.sd.inv = drop(grid.A %*% sd.inv[(1:n.s)+(n.s*t)])
grid.sd.inv = matrix(grid.sd.inv, grid.nx)
image(grid.sd.inv, asp = grid.ratio, col = viridis::inferno(100), zlim = range(sd.inv))

# create gifs
dirname = "img"
width = 720
height = 480
plot.gif(mu.inla, bound, mesh.s, mesh.t, viridis::viridis(200), width, height, dirname, "mu_inla.gif")
plot.gif(sd.inla, bound, mesh.s, mesh.t, viridis::inferno(200), width, height, dirname, "sd_inla.gif")
plot.gif(mu.inv, bound, mesh.s, mesh.t, viridis::viridis(200), width, height, dirname, "mu_inv.gif")
plot.gif(sd.inv, bound, mesh.s, mesh.t, viridis::inferno(200), width, height, dirname, "sd_inv.gif")

