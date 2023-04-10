# command line arguments
args = commandArgs(trailingOnly=TRUE)
for(i in seq_along(args)){
    # problem dimensions
    if(args[i] == "-ms") m.s = as.integer(args[i+1])
    if(args[i] == "-mt") m.t = as.integer(args[i+1])
    if(args[i] == "-res") res.s = as.integer(args[i+1])
    if(args[i] == "-nh") n.h = as.integer(args[i+1])
    
    # hyperparameters
    if(args[i] == "-fixed") fixed = as.integer(args[i+1])
    if(args[i] == "-rt") range.t = as.double(args[i+1])
    if(args[i] == "-rs") range.s = as.double(args[i+1])
    if(args[i] == "-sigma") sigma = as.double(args[i+1])
    if(args[i] == "-tauy") tau.y = as.double(args[i+1])
    if(args[i] == "-taub") tau.b = as.double(args[i+1])
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

# prepare data
source("prepdata.R")

# ------------------------------------------------------------------------------
# ---------------------------- SET HYPERPARAMETERS -----------------------------
# ------------------------------------------------------------------------------

# set smoothness (critical diffusion, 121)
d = 2
alpha.t = 1 # temporal order
alpha.s = 2
alpha.e = d/2
alpha = alpha.e + alpha.s * (alpha.t-1/2) # spatial order
nu.t = alpha.t - 1/2
nu.s = alpha.s * nu.t

# set hyperparameters
if(!exists(deparse(substitute(r.s)))) r.s = 1259
if(!exists(deparse(substitute(r.t)))) r.t = 50
if(!exists(deparse(substitute(sigma.st)))) sigma.st = 5.96
if(!exists(deparse(substitute(tau.y)))) tau.y = 0.143
if(!exists(deparse(substitute(tau_b)))) tau_b = 1e-5
sigma.sq = sigma.st^2

# convert hyperparameters
c.1 = gamma(alpha.t-1/2) / gamma(alpha.t) / (4*pi)^(1/2)
c.2 = gamma(alpha-d/2) / gamma(alpha) / (4*pi)^(d/2)
gamma.s = sqrt(8*nu.s) / range.s
gamma.t = range.t * gamma.s^alpha.s / sqrt(8*(alpha.t-1/2))
gamma.e = sqrt(c.1*c.2 / gamma.t / gamma.s^(2*alpha-d) / sigma.sq)

# ------------------------------------------------------------------------------
# ---------------------------- CREATE MATRICES ---------------------------------
# ------------------------------------------------------------------------------

# data matrix
A.b = as(as.matrix(data[-(1:3)]) + .Machine$double.xmin, "RsparseMatrix")
id.na = which(is.na(data$y))
n.na = length(id.na)

# fem matrices
fem.t = inla.mesh.fem(mesh.t, order = 2)
fem.s = inla.mesh.fem(mesh.s, order = 3)

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

# ------------------------------------------------------------------------------
# ---------------------------- SAVE MATRICES -----------------------------------
# ------------------------------------------------------------------------------

# save matrices
path = "../../data/"
write_petsc_mat(J.0*gamma.e^2, paste0(path,"J0"))
write_petsc_mat(J.1*gamma.e^2*2*gamma.t, paste0(path,"J1"))
write_petsc_mat(J.2*gamma.e^2*gamma.t^2, paste0(path,"J2"))
write_petsc_mat(K.1, paste0(path,"K1"))
write_petsc_mat(K.2, paste0(path,"K2"))
write_petsc_mat(K.3, paste0(path,"K3"))
write_petsc_mat(A.t*sqrt(tau.y), paste0(path,"At"))
write_petsc_mat(A.s, paste0(path,"As"))
write_petsc_mat(A.b*sqrt(tau.y), paste0(path,"Ab"))
write_petsc_vec(y*sqrt(tau.y), paste0(path,"y"))
write_petsc_is(id.na, paste0(path,"isna"))
