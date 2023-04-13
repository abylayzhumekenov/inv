# mu = readBin("../../data/mu", "double", 1+n.st+n.b, endian="swap")[-1]
# d = readBin("../../data/d", "double", 1+n.st+n.b, endian="swap")[-1]

n.long = 1000
n.lat = 600
lattice.long = seq(-11000, -5000, length = n.long)
lattice.lat = seq(2600, 6000, length = n.lat)
lattice.ratio = diff(range(lattice.lat)) / diff(range(lattice.long))
loc.lattice = expand.grid(lattice.long, lattice.lat)
colnames(loc.lattice) = colnames(loc)
loc.lattice = as.matrix(loc.lattice)
A.lattice = inla.spde.make.A(mesh = mesh.s, loc = loc.lattice)

mu.lattice = drop(A.lattice %*% mu[(1:n.s)+(n.s*364)])
mu.lattice = matrix(mu.lattice, n.long)
image(mu.lattice, asp = lattice.ratio, col = viridis::viridis(100))

d.lattice = drop(A.lattice %*% d[(1:n.s)+(n.s*1000)])
d.lattice = matrix(d.lattice, n.long)
image(d.lattice, asp = lattice.ratio, col = rev(viridis::inferno(100)))
