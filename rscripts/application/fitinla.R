library(INLA)
library(INLAspacetime)
library(inlabru)
inla.setOption(smtp = "pardiso", inla.mode = "compact", pardiso.license = "~/pardiso.license")
# inla.setOption(smtp = "taucs", inla.mode = "compact")

# load the data
source("getdata.R")
detach("package:data.table", unload = TRUE)
wdat = wdat[1:10,1:10]

# create spatial and temporal mesh
resolution.s = 1000
resolution.t = 9
bound = inla.nonconvex.hull(points = coordinates(stations), 
                            convex = 200, concave = 200, resolution = 100)
mesh.s = inla.mesh.2d(boundary = bound, 
                      max.edge = c(1, 2) * resolution.s, 
                      offset = c(1e-9, resolution.s * 2), 
                      cutoff = resolution.s / 4)
mesh.t = inla.mesh.1d(1:resolution.t)

# get dimensions
wdat = wdat[,1:(resolution.t+1)]
n.t = mesh.t$n
n.s = mesh.s$n
m.t = dim(wdat)[2]-1
m.s = dim(wdat)[1]
n.st = n.s * n.t
m.st = m.s * m.t

# prepare data
loc = stations@coords[match(wdat$station, stations$station),]
elevation = stations$elevation[match(wdat$station, stations$station)]
data = data.frame(longitude = rep(loc[,1], m.t),
                  latitude = rep(loc[,2], m.t),
                  time = rep(1:m.t, each = m.s),
                  y = c(as.matrix(wdat[,-1])) / 10,
                  elevation = rep(elevation, m.t),
                  # northing = rep(loc[,2], m.t),
                  sine = rep(sin(2*pi*(1:m.t)/365.25), each = m.s),
                  cosine = rep(cos(2*pi*(1:m.t)/365.25), each = m.s))

# define a model
model = ~ -1 + Interecept(1) + 
    # elevation + northing +
    field(list(space = cbind(longitude, latitude),
               time = time),
          model = stmodel)
stmodel = stModel.define(mesh.s, mesh.t, "121",
                         control.priors = list(prs = c(9.895, 0.0),
                                               prt = c(4.026, 0.0),
                                               psigma = c(5.596, 0.0)))
# lkprec = list(prec = list(initial = -1.289, fixed = TRUE))
lkprec = list(prec = list(initial = 10, fixed = TRUE))

# fit the model
res = bru(model,
          like(formula = y ~ .,
               family = "gaussian",
               control.family = list(hyper = lkprec),
               data = data),
          options = list(verbose = TRUE,
                         safe = FALSE,
                         num.threads = "1:1",
                         control.inla = list(int.strategy = "eb"),
                         control.compute = list(config = TRUE)))

print(res$summary.hyperpar)
