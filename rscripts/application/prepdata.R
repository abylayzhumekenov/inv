# get the data
source("getdata.R")
detach("package:data.table", unload = TRUE)
set.seed(1)

# set problem dimensions
if(!exists(deparse(substitute(m.t)))) m.t = dim(wdat)[2]-1
if(!exists(deparse(substitute(m.s)))) m.s = dim(wdat)[1]
if(!exists(deparse(substitute(res.s)))) res.s = 500
if(!exists(deparse(substitute(n.h)))) n.h = 1

# create mesh
mesh.t = inla.mesh.1d(1:m.t)
bound = inla.nonconvex.hull(points = stations@coords, convex = 200, concave = 200, resolution = 100)
mesh.s = inla.mesh.2d(boundary = bound, max.edge = c(1, 2)*res.s, offset = c(1e-3, 3*res.s), cutoff = res.s/4)
n.t = mesh.t$n
n.s = mesh.s$n

# prepare data
wdat = wdat[sort(sample(1:(dim(wdat)[1]), m.s)), 1:(m.t+1)]
loc = stations@coords[match(wdat$station, stations$station),]
ele = stations$elevation[match(wdat$station, stations$station)]
data = data.frame(xcoord = rep(loc[,1], m.t),
                  ycoord = rep(loc[,2], m.t),
                  time = rep(1:m.t, each = m.s),
                  y = c(as.matrix(wdat[,-1])) / 10,
                  elevation = rep(ele, m.t) / 1000,
                  northing = rep(loc[,2], m.t) / 10000)
for(i in 1:n.h){
    harm = data.frame(sin = rep(sin(i*2*pi*(1:m.t-1)/365.25), each = m.s),
                      cos = rep(cos(i*2*pi*(1:m.t-1)/365.25), each = m.s))
    names(harm) = c(paste0("harmonic", i, ".sin"), paste0("harmonic", i, ".cos"))
    data = cbind(data, harm)
}





