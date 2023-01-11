# command line arguments
args = commandArgs(trailingOnly=TRUE)
for(i in seq_along(args)){
	if(args[i] == "-mt") m.t = as.integer(args[i+1])
	if(args[i] == "-ms") m.s = as.integer(args[i+1])
	if(args[i] == "-rt") range.t = as.integer(args[i+1])
}

# load libraries
library(INLA)
set.seed(100)

# function to write in binary
write_petsc = function(Q, filename){
    # encode the matrix
    Q = as(Q, "generalMatrix")
    x = list(classid = 1211216,
             nrows = nrow(Q),
             ncols = nrow(Q),
             nnz = length(Q@x),
             nnz_row = diff(Q@p),
             nnz_i = Q@i,
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
sigma.sq = 1
range.s = 1
if(!exists(deparse(substitute(range.t)))) range.t = 10
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
if(!exists(deparse(substitute(m.t)))) m.t = 100
mesh.t = inla.mesh.1d(1:m.t)
fem.t = inla.mesh.fem(mesh.t, order = 2)

# spatial mesh
if(!exists(deparse(substitute(m.s)))) m.s = 4
# loc.s = cbind(rep(seq(0, 1, length.out = m.s), times = m.s), rep(seq(0, 1, length.out = m.s), each = m.s))
# mesh.s = inla.mesh.2d(loc.s, max.edge = 2/m.s)
mesh.s = inla.mesh.create(globe = m.s)
fem.s = inla.mesh.fem(mesh.s, order = 3)

# temporal matrices
J.0 = fem.t$c0
J.1 = Diagonal(m.t, c(0.5, rep(0, m.t-2), 0.5))
J.2 = fem.t$g1

# spatial matrices
K.1 = gamma.s^2 * fem.s$c0 + fem.s$g1
K.2 = gamma.s^4 * fem.s$c0 + 2*gamma.s^2 * fem.s$g1 + fem.s$g2
K.3 = gamma.s^6 * fem.s$c0 + 3*gamma.s^4 * fem.s$g1 + 3*gamma.s^2 * fem.s$g2 + fem.s$g3

# save matrices
write_petsc(J.0*gamma.e^2, "J0")
write_petsc(J.1*gamma.e^2*2*gamma.t, "J1")
write_petsc(J.2*gamma.e^2*gamma.t^2, "J2")
write_petsc(K.1, "K1")
write_petsc(K.2, "K2")
write_petsc(K.3, "K3")

print(c(m.t, m.s, range.t))
# # spatio-temporal precision
# Q.st = gamma.e^2 * (kronecker(J.0, K.3) + kronecker(J.1*2*gamma.t, K.2) + kronecker(J.2*gamma.t^2, K.1))
# n.st = nrow(Q.st)
# n.t = m.t
# n.s = n.st / n.t
# 
# # fixed effects precision
# n.beta = 1
# tau.beta = 1e-5
# Q.beta = Diagonal(n.beta, tau.beta)
# 
# # latent field prior precision
# Q.prior = rbind(cbind(Q.st, sparseMatrix(NULL, NULL, dims = c(nrow(Q.st), nrow(Q.beta)))),
#                 cbind(sparseMatrix(NULL, NULL, dims = c(nrow(Q.beta), nrow(Q.st))), Q.beta))
# n.x = nrow(Q.prior)
# 
# # data matrices
# Z = matrix(rnorm(n.st*n.beta), n.st, n.beta)
# B = Diagonal(n.st)
# A = cbind(B, Z)
# 
# # observation precision
# tau.y = 1e-5
# Q.y = Diagonal(n.st, tau.y)
# 
# # latent field posterior precision
# Q.posterior = Q.prior + crossprod(A, Q.y) %*% A
# 
# # # ------------------------------------------------------------------------------
# # # PLOT PRECISION MATRICES
# # # ------------------------------------------------------------------------------
# # 
# # # plot the prior precision
# # pdf("Q_prior.pdf")
# # image(Q.prior!=0, border.col = NA)
# # dev.off()
# # 
# # # plot the posterior precision
# # pdf("Q_posterior.pdf")
# # image(Q.posterior!=0, border.col = NA)
# # dev.off()
# 
# # ------------------------------------------------------------------------------
# # INVERSION USING R
# # ------------------------------------------------------------------------------
# 
# # compute inverse
# r.time = system.time(d.true <- diag(solve(Q.posterior[1:n.st, 1:n.st])))
# cat(paste("Time elapsed:", formatC(r.time[[3]], format = "e"), "sec"), "\n")
# 
# # ------------------------------------------------------------------------------
# # CALL C CODE FOR INVERSION
# # ------------------------------------------------------------------------------
# 
# # set C options
# c.options = list(n_core = 2,
#                  n_sample = 100,
#                  n_niter = 1000,
#                  n_neighbor = 10,
#                  verbose = 1,
#                  verbose_s = 0)
# 
# # # save the posterior precision
# # writeMM(Q.posterior[1:n.st, 1:n.st], c.options$fin)
# 
# # call C code
# c.call = paste(paste0("PETSC_DIR=", Sys.getenv("PETSC_DIR")),
#                paste0("PETSC_ARCH=", Sys.getenv("PETSC_ARCH")),
#                paste0("METIS_DIR=", Sys.getenv("METIS_DIR")),
#                paste0("LD_LIBRARY_PATH=", Sys.getenv("LD_LIBRARY_PATH")),
#                "export PETSC_DIR",
#                "export PETSC_ARCH",
#                "export METIS_DIR",
#                "export LD_LIBRARY_PATH",
#                "cd ../",
#                paste("mpiexec",
#                      "-n", c.options$n_core,
#                      "./bin/main",
#                      "-ns", c.options$n_sample,
#                      "-nmax", c.options$n_niter,
#                      "-nn", c.options$n_neighbor,
#                      "-v", c.options$verbose,
#                      "-vs", c.options$verbose_s),
#                sep = "\n")
# c.time = system.time(c.out <- system(command = c.call, intern = TRUE, wait = TRUE))
# cat(paste(c.out, collapse = "\n"), "\n")
# cat(paste("Time elapsed:", formatC(c.time[[3]], format = "e"), "sec"), "\n")
# 
# # read the result
# d.inv = readBin("out", "double", n.st + 1, endian = "swap")[-1]
# 
# # ------------------------------------------------------------------------------
# # COMPARE
# # ------------------------------------------------------------------------------
# 
# # compare results
# rmse = abs(d.inv/d.true-1)
# plot(rmse*0, t = "l", ylim = c(0,1))
# lines(rmse, col = "red")
# legend("topright", legend = c("Krylov", "True"), col = c("red", "black"), lty = 1)
# d.error = norm(d.inv - d.true, "2")
# 
# # # ------------------------------------------------------------------------------
# # # OBTAIN POSTERIOR USING R
# # # ------------------------------------------------------------------------------
# 
# # # generate observations
# # x = drop(inla.qsample(1, Q.prior))
# # y = drop(A %*% x + inla.qsample(1, Q.y))
# 
# # # compute the posterior mean
# # mu.x = drop(solve(Q.posterior, crossprod(A, Q.y %*% y)))
# # mu.s = mu.x[1:n.s]
# 
# # # ------------------------------------------------------------------------------
# # # PROJECT SOLUTION AND PLOT
# # # ------------------------------------------------------------------------------
# 
# # # project to sphere
# # long.pred = seq(-180, 180, length.out = 360*1)
# # lat.pred = seq(-90, 90, length.out = 180*1)
# # loc.pred = as.matrix(expand.grid(long.pred, lat.pred))
# # loc.pred.sphere = inla.mesh.map(loc.pred, "longlat")
# # A.pred = inla.spde.make.A(mesh.s, loc.pred.sphere)
# 
# # # compute interpolated mean
# # mu.pred = drop(A.pred %*% mu.s)
# # mu.pred.matrix = matrix(mu.pred, nrow = length(lat.pred), byrow = TRUE)
# 
# # # # visualize the posterior mean
# # # image(t(mu.pred.matrix), ylim = c(1,0), asp = 0.5, col = hcl.colors(100, "Inferno"))
# 
# # ------------------------------------------------------------------------------
# # PRINT RESULTS
# # ------------------------------------------------------------------------------
# 
# cat(paste0("Dimensions|", "\t\tSpace: ", n.s, "\tTime: ", n.t, "\tOverall: ", n.st, "\n",
#            "Number of cores|", "\tR code: ", 1, "\t\tC code: ", c.options$n_core, "\n",
#            "Time elapsed|", "\t\tR code: ", formatC(r.time[[3]], format = "e"),
#            "\tC code: ", formatC(c.time[[3]], format = "e")), "\n")
