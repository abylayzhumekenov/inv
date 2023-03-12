library(INLA)
library(INLAspacetime)
library(inlabru)
inla.setOption(smtp='taucs', inla.mode='compact')

# load the data
source("getdata.R")
detach("package:data.table", unload = TRUE)
set.seed(1)
# wdat = wdat[sample(dim(wdat)[1], 100), 1:20]
wdat = wdat[, 1:(365*3)]

# temporal and spatial mesh
mesh.t = inla.mesh.1d(1:(dim(wdat)[2]-1))
loc = stations@coords[match(wdat$station, stations$station),]
elevation = stations$elevation[match(wdat$station, stations$station)]
bound = inla.nonconvex.hull(points = loc, convex = 200, concave = 200, resolution = 100)
mesh.s = inla.mesh.create(loc = loc, boundary = bound, cutoff = 200)

# set dimensions
n.t = mesh.t$n
n.s = mesh.s$n
m.t = dim(wdat)[2]-1
m.s = dim(wdat)[1]
n.st = n.s * n.t
m.st = m.s * m.t

# prepare data
data = data.frame(longitude = rep(loc[,1], m.t),
                  latitude = rep(loc[,2], m.t),
                  time = rep(1:m.t, each = m.s),
                  y = c(as.matrix(wdat[,-1])),
                  elevation = rep(elevation, m.t),
                  northing = rep(loc[,2], m.t),
                  sine = rep(sin(2*pi*(1:m.t)/365.25), each = m.s),
                  cosine = rep(cos(2*pi*(1:m.t)/365.25), each = m.s))

# define a model
model = ~ -1 + Interecept(1) + elevation + northing + 
    field(list(space = cbind(longitude, latitude),
               time = time),
          model = stmodel)
stmodel = stModel.define(mesh.s, mesh.t, "121",
                         control.priors = list(prs = c(100, 0.01),
                                               prt = c(1, 0.01),
                                               psigma = c(10000, 0.01)))
lkprec = list(prec = list(initial = 10, fixed = FALSE))

# fit the model
res = bru(model,
          like(formula = y ~ .,
               family = "gaussian",
               control.family = list(hyper = lkprec),
               data = data),
          options = list(verbose = TRUE,
                         num.threads = "5:6",
                         control.inla = list(int.strategy = "eb"),
                         control.compute = list(config = TRUE)))

print(res$summary.hyperpar)
