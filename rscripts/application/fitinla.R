library(INLA)
library(INLAspacetime)
library(inlabru)
inla.setOption(smtp = "pardiso", inla.mode = "compact", pardiso.license = "~/pardiso.license")

# load the data
source("getdata.R")
detach("package:data.table", unload = TRUE)

# set dimensions
m.t = n.t = 365*4+1
m.s = 1000
m.st = m.t * m.s

# create mesh
mesh.t = inla.mesh.1d(1:n.t)
resolution.s = 200
bound = inla.nonconvex.hull(points = stations@coords, convex = 200, concave = 200, resolution = 100)
mesh.s = inla.mesh.2d(boundary = bound, max.edge = c(1, 2)*resolution.s, 
                      offset = c(1e-3, 3*resolution.s), cutoff = resolution.s/4)
n.s = mesh.s$n
n.st = n.t * n.s

# prepare data
wdat = wdat[sort(sample(1:(dim(wdat)[1]), m.s)), 1:(m.t+1)]
loc = stations@coords[match(wdat$station, stations$station),]
ele = stations$elevation[match(wdat$station, stations$station)]
data = data.frame(longitude = rep(loc[,1], m.t),
                  latitude = rep(loc[,2], m.t),
                  time = rep(1:m.t, each = m.s),
                  y = c(as.matrix(wdat[,-1])) / 10,
                  elevation = rep(ele, m.t) / 1000,
                  northing = rep(loc[,2], m.t) / 10000,
                  harmonic1 = rep(sin(2*pi*(1:m.t-1)/365.25), each = m.s),
                  harmonic2 = rep(cos(2*pi*(1:m.t-1)/365.25), each = m.s))

# define a model
model = ~ -1 + Intercept(1) + elevation + northing + harmonic1 + harmonic2 +
    field(list(space = cbind(latitude, longitude), time = time), 
          model = stmodel)
# theta.hat = c(-1.289, 9.895, 14.026, 5.596)
stmodel = stModel.define(mesh.s, mesh.t, "121", 
                         control.priors = list(prs = c(7, 0.1),
                                               prt = c(3, 0.1),
                                               psigma = c(3, 0.1)))
lkprec = list(prec = list(initial = 10, fixed = FALSE))

# fit the model
result = bru(model, 
             like(formula = y ~ ., 
                  family = "gaussian",
                  control.family = list(hyper = lkprec), 
                  data = data),
             options = list(verbose = TRUE,
                            safe = FALSE,
                            num.threads = "8:4",
                            control.inla = list(int.strategy = "eb")))

# print the result
print(result$summary.fixed)
print(result$summary.hyperpar)
