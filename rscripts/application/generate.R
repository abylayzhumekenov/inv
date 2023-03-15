# command line arguments
args = commandArgs(trailingOnly=TRUE)
for(i in seq_along(args)){
	if(args[i] == "-rt") range.t = as.integer(args[i+1])
	if(args[i] == "-rs") range.t = as.integer(args[i+1])
	if(args[i] == "-ss") sigma.sq = as.integer(args[i+1])
	if(args[i] == "-tauy") tau_y = as.double(args[i+1])
	if(args[i] == "-taub") tau_b = as.double(args[i+1])
}

# function to write in binary
write_petsc_mat = function(Q, filename){
    # encode the matrix
    Q = as(Q, "RsparseMatrix")
    x = list(classid = 1211216,
             nrows = nrow(Q),
             ncols = ncol(Q),
             nnz = length(Q@x),
             nnz_row = diff(Q@p),
             nnz_i = Q@j,
             nnz_val = Q@x)
    
    # write the binary data
    fwrite = file(filename, "wb")
    writeBin(object = as.integer(x$classid), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nrows), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$ncols), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nnz), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nnz_row), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nnz_i), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.double(x$nnz_val), con = fwrite, size = 8, endian = "swap")
    close(con = fwrite)
}

write_petsc_vec = function(y, filename){
    # encode the vector
    x = list(classid = 1211214,
             nrows = length(y),
             val = drop(y))
    
    # write the binary data
    fwrite = file(filename, "wb")
    writeBin(object = as.integer(x$classid), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nrows), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.double(x$val), con = fwrite, size = 8, endian = "swap")
    close(con = fwrite)
}

write_petsc_is = function(is, filename){
    # encode the index set
    x = list(classid = 1211218,
             nrows = length(is),
             val = drop(is)-1)
    
    # write the binary data
    fwrite = file(filename, "wb")
    writeBin(object = as.integer(x$classid), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nrows), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$val), con = fwrite, size = 4, endian = "swap")
    close(con = fwrite)
}

# load libraries
library(INLA)
set.seed(100)
source("getdata.R")

# ------------------------------------------------------------------------------
# SET MODEL HYPERPARAMETERS
# ------------------------------------------------------------------------------

# set smoothness (critical diffusion)
d = 2
alpha.t = 1 # temporal order
alpha.s = 2
alpha.e = d/2
alpha = alpha.e + alpha.s * (alpha.t-1/2) # spatial order
nu.t = alpha.t - 1/2
nu.s = alpha.s * nu.t

# set practical range and marginal variance
if(!exists(deparse(substitute(sigma.sq)))) sigma.sq = 4054.745
if(!exists(deparse(substitute(range.s)))) range.s = 1500.000
if(!exists(deparse(substitute(range.t)))) range.t = 2
theta = log(c(sigma.sq, range.s, range.t))

# convert hyperparameters
c.1 = gamma(alpha.t-1/2) / gamma(alpha.t) / (4*pi)^(1/2)
c.2 = gamma(alpha-d/2) / gamma(alpha) / (4*pi)^(d/2)
gamma.s = sqrt(8*nu.s) / range.s
gamma.t = range.t * gamma.s^alpha.s / sqrt(8*(alpha.t-1/2))
gamma.e = sqrt(c.1*c.2 / gamma.t / gamma.s^(2*alpha-d) / sigma.sq)

# ------------------------------------------------------------------------------
# CREATE MESH, MATRICES AND DATA
# ------------------------------------------------------------------------------

# temporal mesh
mesh.t = inla.mesh.1d(1:(dim(wdat)[2]-1))
fem.t = inla.mesh.fem(mesh.t, order = 2)

# spatial mesh
loc = stations@coords[match(wdat$station, stations$station),]
ele = stations$elevation[match(wdat$station, stations$station)]
bound = inla.nonconvex.hull(points = loc, convex = 200, concave = 200, resolution = 100)
mesh.s = inla.mesh.create(loc = loc, boundary = bound)
# mesh.s = inla.mesh.2d(loc = loc, boundary = bound, max.edge = 200)
fem.s = inla.mesh.fem(mesh.s, order = 3)

# set dimensions
n.t = dim(fem.t$g1)[1]
n.s = dim(fem.s$g1)[1]
m.t = dim(wdat)[2] - 1
m.s = dim(wdat)[1]
n.st = n.t * n.s
m.st = m.t * m.s
n.na = sum(is.na(wdat[,-1]))

# set precisions
if(!exists(deparse(substitute(tau_y)))) tau_y = 1e-2
if(!exists(deparse(substitute(tau_b)))) tau_b = 1e-5

# temporal matrices
J.0 = fem.t$c0
J.1 = Diagonal(n.t, c(0.5, rep(0, n.t-2), 0.5))
J.2 = fem.t$g1

# spatial matrices
K.1 = gamma.s^2 * fem.s$c0 + fem.s$g1
K.2 = gamma.s^4 * fem.s$c0 + 2*gamma.s^2 * fem.s$g1 + fem.s$g2
K.3 = gamma.s^6 * fem.s$c0 + 3*gamma.s^4 * fem.s$g1 + 3*gamma.s^2 * fem.s$g2 + fem.s$g3

# projection matrices
A.t = inla.spde.make.A(mesh = mesh.t)
A.s = inla.spde.make.A(mesh = mesh.s, loc = loc)

# observations and covariates
n.h = 4
y = c(as.matrix(wdat[,-1])) / 10
intercept = rep(1, m.st)
elevation = rep(ele, m.t) / 1000
northing = rep(loc[,2], m.t) / 10000
A.b = cbind(intercept, elevation, northing)
for(i in 1:n.h){
    A.b = cbind(A.b, rep(sin(i*2*pi*(1:m.t-1)/365.25), each = m.s))
    A.b = cbind(A.b, rep(cos(i*2*pi*(1:m.t-1)/365.25), each = m.s))
}
A.b = as(A.b + .Machine$double.xmin, "RsparseMatrix")
id.na = which(is.na(y))

# save matrices
path = "../../data/"
write_petsc_mat(J.0*gamma.e^2, paste0(path,"J0"))
write_petsc_mat(J.1*gamma.e^2*2*gamma.t, paste0(path,"J1"))
write_petsc_mat(J.2*gamma.e^2*gamma.t^2, paste0(path,"J2"))
write_petsc_mat(K.1, paste0(path,"K1"))
write_petsc_mat(K.2, paste0(path,"K2"))
write_petsc_mat(K.3, paste0(path,"K3"))
write_petsc_mat(A.t*sqrt(tau_y), paste0(path,"At"))
write_petsc_mat(A.s, paste0(path,"As"))
write_petsc_mat(A.b*sqrt(tau_y), paste0(path,"Ab"))
write_petsc_vec(y*sqrt(tau_y), paste0(path,"y"))
write_petsc_is(id.na, paste0(path,"isna"))
