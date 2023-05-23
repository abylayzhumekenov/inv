# ------------------------------------------------------------------------------
# ---------- SET PARAMETERS ----------------------------------------------------
# ------------------------------------------------------------------------------

# command line arguments
args = commandArgs(trailingOnly=TRUE)
for(i in seq_along(args)){
    if(args[i] == "-ns") n.s = as.integer(args[i+1])
    if(args[i] == "-nt") n.t = as.integer(args[i+1])
    if(args[i] == "-ms") m.s = as.integer(args[i+1])
    if(args[i] == "-mt") m.t = as.integer(args[i+1])
    if(args[i] == "-res1") res1 = as.double(args[i+1])
    if(args[i] == "-res2") res2 = as.double(args[i+1])
    if(args[i] == "-res3") res3 = as.double(args[i+1])
    if(args[i] == "-rs") r.s = as.double(args[i+1])
    if(args[i] == "-rt") r.t = as.double(args[i+1])
    if(args[i] == "-sig") sig = as.double(args[i+1])
    if(args[i] == "-tauy") tau.y = as.double(args[i+1])
    if(args[i] == "-taub") tau.b = as.double(args[i+1])
}

# default hyperparameter values
if(!exists(deparse(substitute(r.s)))) r.s = 1259
if(!exists(deparse(substitute(r.t)))) r.t = 50
if(!exists(deparse(substitute(sig)))) sig = 5.96
if(!exists(deparse(substitute(tau.y)))) tau.y = 0.143
if(!exists(deparse(substitute(tau.b)))) tau.b = 1e-3

# ------------------------------------------------------------------------------
# ---------- PREPARE THE DATA --------------------------------------------------
# ------------------------------------------------------------------------------

# load libraries
suppressMessages(suppressWarnings(library(INLA)))
source("getdata.R")
detach("package:data.table", unload = TRUE)
rm(list = setdiff(ls(), c("wdat", "stations", "n.s", "n.t", "m.s", "m.t",
                          "res1", "res2", "res3", "r.s", "r.t", "sig", "tau.y", "tau.b")))
source("write_petsc.R")
set.seed(1)

# set mesh resolution and data dimensions
if(!exists(deparse(substitute(res1)))) res1 = 200
if(!exists(deparse(substitute(res2)))) res2 = 500
if(!exists(deparse(substitute(res3)))) res3 = 2000
if(!exists(deparse(substitute(m.s)))) m.s = dim(wdat)[1]
if(!exists(deparse(substitute(m.t)))) m.t = dim(wdat)[2]-1

# create mesh
bound = inla.nonconvex.hull(points = stations@coords, convex = res1, concave = res1, resolution = 100)
mesh.s = inla.mesh.2d(boundary = bound, max.edge = c(1, 4)*res2, offset = c(1e-3, res3), cutoff = res2/4)
mesh.t = inla.mesh.1d(1:m.t)
n.s = mesh.s$n
n.t = mesh.t$n

# covariates
wdat = wdat[sort(sample(1:(dim(wdat)[1]), m.s)), 1:(m.t+1)]
loc = stations@coords[match(wdat$station, stations$station),]
ele = stations$elevation[match(wdat$station, stations$station)]
A.b = cbind(intercept = rep(1, m.s*m.t),
            elevation = rep(ele, m.t) / 1000,
            northing = rep(loc[,2], m.t) / 10000,
            sin = rep(sin(2*pi*(1:m.t-1)/365.25), each = m.s),
            cos = rep(cos(2*pi*(1:m.t-1)/365.25), each = m.s))
n.b = ncol(A.b)

# observations
y = c(as.matrix(wdat[,-1])) / 10
y.na = is.na(y)
y[y.na] = 0

# ------------------------------------------------------------------------------
# ---------- SET HYPERPARAMETERS -----------------------------------------------
# ------------------------------------------------------------------------------

# set smoothness (121, critical diffusion)
d.s = 2
alpha.t = 1 # temporal order
alpha.s = 2
alpha.e = d.s/2
alpha = alpha.e + alpha.s * (alpha.t-1/2) # spatial order
nu.t = alpha.t - 1/2
nu.s = alpha.s * nu.t

# convert hyperparameters
c.1 = gamma(alpha.t-1/2) / gamma(alpha.t) / (4*pi)^(1/2)
c.2 = gamma(alpha-d.s/2) / gamma(alpha) / (4*pi)^(d.s/2)
gamma.s = sqrt(8*nu.s) / r.s
gamma.t = r.t * gamma.s^alpha.s / sqrt(8*(alpha.t-1/2))
gamma.e = sqrt(c.1*c.2 / gamma.t / gamma.s^(2*alpha-d.s) / sig^2)

# ------------------------------------------------------------------------------
# ---------- GENERATE PRECISION MATRICES ---------------------------------------
# ------------------------------------------------------------------------------

# create fem objects
fem.t = inla.mesh.fem(mesh.t, order = 2)
fem.s = inla.mesh.fem(mesh.s, order = 3)

# temporal matrices
J.0 = gamma.e^2 * fem.t$c0
J.1 = gamma.e^2 * gamma.t * Diagonal(n.t, c(0.5, rep(0, n.t-2), 0.5))
J.2 = gamma.e^2 * gamma.t^2 * fem.t$g1

# spatial matrices
K.1 = gamma.s^2 * fem.s$c0 + fem.s$g1
K.2 = gamma.s^4 * fem.s$c0 + 2*gamma.s^2 * fem.s$g1 + fem.s$g2
K.3 = gamma.s^6 * fem.s$c0 + 3*gamma.s^4 * fem.s$g1 + 3*gamma.s^2 * fem.s$g2 + fem.s$g3

# projection matrices
A.t = inla.spde.make.A(mesh = mesh.t)
A.s = inla.spde.make.A(mesh = mesh.s, loc = loc)

# observation precision and covariates
q.yy = sqrt(tau.y) * (!y.na) # sqrt(Q.y)

# ------------------------------------------------------------------------------
# ---------- SAVE THE OBJECTS --------------------------------------------------
# ------------------------------------------------------------------------------

# save matrices
path = "../../data/"
write_petsc_mat(J.0, paste0(path, "J0"))
write_petsc_mat(J.1, paste0(path, "J1"))
write_petsc_mat(J.2, paste0(path, "J2"))
write_petsc_mat(K.1, paste0(path, "K1"))
write_petsc_mat(K.2, paste0(path, "K2"))
write_petsc_mat(K.3, paste0(path, "K3"))
write_petsc_mat(t(A.t), paste0(path, "At"))
write_petsc_mat(t(A.s), paste0(path, "As"))
write_petsc_mat(A.b + .Machine$double.xmin, paste0(path, "Ab"))
write_petsc_vec(y, paste0(path, "y"))
write_petsc_vec(q.yy, paste0(path, "qyy"))

# save hyperparameters and mesh info for later use
save(list = c("n.s", "n.t", "m.s", "m.t", "n.b",
              "r.s", "r.t", "sig", "tau.y", "tau.b",
              "mesh.s", "mesh.t", "bound", "y", "loc", "ele"), file = "data/data.Rdata")

# print out info
cat(format(c("ns","nt","nu","ms","mt","mu","nb"), width=12, justify="right"), "\n",
    format(c(n.s,n.t,n.s*n.t,m.s,m.t,m.s*m.t,n.b), width=12, justify="right"), "\n",
    "\n",
    format(c("rs","rt","sig","tauy","taub"), width=12, justify="right"), "\n",
    format(c(r.s,r.t,sig,tau.y,tau.b), width=12, justify="right"), "\n",
    sep = "")
