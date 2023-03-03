# command line arguments
args = commandArgs(trailingOnly=TRUE)
for(i in seq_along(args)){
	if(args[i] == "-mt") m.t = as.integer(args[i+1])
	if(args[i] == "-ms") m.s = as.integer(args[i+1])
	if(args[i] == "-rt") range.t = as.integer(args[i+1])
	if(args[i] == "-rs") range.t = as.integer(args[i+1])
	if(args[i] == "-ss") sigma.sq = as.integer(args[i+1])
}

# generate data
source("generate.R")

# if the problem is too large, abort
if(n.st > 1000) {
    cat("Problem is too large, try smaller than 1000\n")
    quit()
}

# prepare matrices
J.0 = J.0*gamma.e^2
J.1 = J.1*gamma.e^2*2*gamma.t
J.2 = J.2*gamma.e^2*gamma.t^2
A.t = A.t*sqrt(tau_y)
A.b = A.b*sqrt(tau_y)
y = y*sqrt(tau_y)

# compute the projection matrix
A.st = kronecker(A.t, A.s)
A = cbind(A.st, A.b)
A[id.na,] = 0
y[id.na] = 0

# compute the posterior precision
Q.posterior = kronecker(J.0, K.3) + kronecker(J.1, K.2) + kronecker(J.2, K.1)
Q.posterior = rbind(cbind(Q.posterior, Matrix(matrix(0, n.st, 1), sparse = TRUE)),
                    cbind(Matrix(matrix(0, 1, n.st), sparse = TRUE), tau_b))
Q.posterior = Q.posterior + t(A) %*% A

# solve using R
mu0 = drop(solve(Q.posterior, t(A) %*% y))
d0 = diag(solve(Q.posterior))

# read PETSc solution
mu = readBin("../../data/mu", "double", n.st+n.b+1, endian = "swap")[-1]
d = readBin("../../data/d", "double", n.st+n.b+1, endian = "swap")[-1]

# print the error norm
cat(paste0("|e_mu|=", norm(mu-mu0, "2"), "\n"))
cat(paste0("|e_d|=", norm(d-d0, "2"), "\n"))

# plot the error graphs
par(mfrow = c(2,1))
plot(mu-mu0, t="l")
abline(h=0, lty=3)
plot(d-d0, t="l")
abline(h=0, lty=3)
par(mfrow = c(1,1))


