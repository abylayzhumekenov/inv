# command line arguments
args = commandArgs(trailingOnly=TRUE)
for(i in seq_along(args)){
	if(args[i] == "-mt") m.t = as.integer(args[i+1])
	if(args[i] == "-ms") m.s = as.integer(args[i+1])
	
	if(args[i] == "-rt") range.t = as.integer(args[i+1])
	if(args[i] == "-rs") range.t = as.integer(args[i+1])
	if(args[i] == "-sigma") sigma.st = as.integer(args[i+1])
	if(args[i] == "-tauy") tau.y = as.double(args[i+1])
	if(args[i] == "-taub") tau.b = as.double(args[i+1])
}

# generate data
source("fitinla.R")

# # if the problem is too large, abort
# if(n.st > 1000) {
#     cat("Problem is too large, try smaller than 1000\n")
#     quit()
# }

# # prepare matrices
# J.0 = J.0*gamma.e^2
# J.1 = J.1*gamma.e^2*2*gamma.t
# J.2 = J.2*gamma.e^2*gamma.t^2
# A.t = A.t*sqrt(tau_y)
# A.b = A.b*sqrt(tau_y)
# y = y*sqrt(tau_y)
# 
# # compute the projection matrix
# A.st = kronecker(A.t, A.s)
# A = cbind(A.st, A.b)
# A[id.na,] = 0
# y[id.na] = 0

# compute the posterior precision
A = cbind(kronecker(A.t, A.s), A.b)
Q.posterior = gamma.e^2 * (kronecker(J.0, K.3) + kronecker(J.1*gamma.t, K.2) + kronecker(J.2*gamma.t^2, K.1))
Q.posterior = rbind(cbind(Q.posterior, Matrix(matrix(0, n.st, n.b), sparse = TRUE)),
                    cbind(Matrix(matrix(0, n.b, n.st), sparse = TRUE), Diagonal(n.b, tau.b)))
Q.posterior = Q.posterior + crossprod(A, Diagonal(n.st, tau.y)) %*% A

# solve using R
mu0 = drop(solve(Q.posterior, crossprod(A, Diagonal(n.st, tau.y)) %*% y))
d0 = diag(solve(Q.posterior))

# solve using INLA
mu1 = c(result$summary.random$field$mode, result$summary.fixed$mode)
d1 = c(result$summary.random$field$sd^2, result$summary.fixed$sd^2)

# read PETSc solution
mu2 = readBin("../../data/mu", "double", n.st+n.b+1, endian = "swap")[-1]
d2 = readBin("../../data/d", "double", n.st+n.b+1, endian = "swap")[-1]

# print the error norm
cat(paste0("|e_mu|=", norm(mu1-mu0, "2"), "\n"))
cat(paste0("|e_d|=", norm(d1-d0, "2"), "\n"))
cat(paste0("|e_mu|=", norm(mu2-mu0, "2"), "\n"))
cat(paste0("|e_d|=", norm(d2-d0, "2"), "\n"))

# plot the error graphs
plot(mu1-mu0, t="l")
lines(mu2-mu0, col="gray")
abline(h=0, lty=3)
plot(d1-d0, t="l")
lines(d2-d0, col="gray")
abline(h=0, lty=3)


