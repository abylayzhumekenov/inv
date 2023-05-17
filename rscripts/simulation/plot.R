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
Q.prior = gamma.e^2 * (kronecker(J.0, K.3) + kronecker(J.1*gamma.t, K.2) + kronecker(J.2*gamma.t^2, K.1))
Q.prior = rbind(cbind(Q.prior, Matrix(matrix(0, n.st, n.b), sparse = TRUE)),
                cbind(Matrix(matrix(0, n.b, n.st), sparse = TRUE), Diagonal(n.b, tau.b)))
Q.posterior = Q.prior + crossprod(A, Diagonal(n.st, tau.y)) %*% A

# solve using R
t0 = Sys.time()
mu0 = drop(solve(Q.posterior, crossprod(A, Diagonal(n.st, tau.y)) %*% y))
d0 = diag(solve(Q.posterior))
sd0 = sqrt(d0)
t0 = Sys.time() - t0

# solve using INLA
mu1 = c(result$summary.random$field$mode, result$summary.fixed$mode)
d1 = c(result$summary.random$field$sd^2, result$summary.fixed$sd^2)
sd1 = sqrt(d1)

# read PETSc solution
mu2 = readBin("../../data/mu", "double", n.st+n.b+1, endian = "swap")[-1]
d2 = readBin("../../data/d", "double", n.st+n.b+1, endian = "swap")[-1]
sd2 = sqrt(d2)

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

# plot projections?..
plot.res = c(128, 128)
coords.spherical = expand.grid(seq(0,2*pi,length=plot.res[1]+1)[-1], seq(0,pi,length=plot.res[2]+1)[-1])
coords.cartesian = data.frame(x = sin(coords.spherical$Var2) * cos(coords.spherical$Var1),
                              y = sin(coords.spherical$Var2) * sin(coords.spherical$Var1),
                              z = cos(coords.spherical$Var2))
A.spherical = inla.spde.make.A(mesh.s, as.matrix(coords.cartesian))

mu0.spherical = drop(A.spherical %*% mu0[1:n.s])
mu0.spherical = matrix(mu0.spherical, plot.res[1])
mu1.spherical = drop(A.spherical %*% mu1[1:n.s])
mu1.spherical = matrix(mu1.spherical, plot.res[1])
mu2.spherical = drop(A.spherical %*% mu2[1:n.s])
mu2.spherical = matrix(mu2.spherical, plot.res[1])

sd0.spherical = drop(A.spherical %*% sd0[1:n.s])
sd0.spherical = matrix(sd0.spherical, plot.res[1])
sd1.spherical = drop(A.spherical %*% sd1[1:n.s])
sd1.spherical = matrix(sd1.spherical, plot.res[1])
sd2.spherical = drop(A.spherical %*% sd2[1:n.s])
sd2.spherical = matrix(sd2.spherical, plot.res[1])

zlim.mu = range(mu0.spherical, mu1.spherical, mu2.spherical)
zlim.sd = range(sd0.spherical, sd1.spherical, sd2.spherical)

pdf("fig.sim.1.pdf", width=8, height=5)
par(mfrow=c(2,3), mar=c(3,3,1,1))
image(mu0.spherical, col=viridis::viridis(100), zlim=zlim.mu, xaxt="n", yaxt="n")
title(ylab="Mean", line=1)
image(mu1.spherical, col=viridis::viridis(100), zlim=zlim.mu, xaxt="n", yaxt="n")
image(mu2.spherical, col=viridis::viridis(100), zlim=zlim.mu, xaxt="n", yaxt="n")
image(sd0.spherical, col=viridis::inferno(100), zlim=zlim.sd, xaxt="n", yaxt="n")
title(xlab="R (base)", ylab="SD", line=1)
image(sd1.spherical, col=viridis::inferno(100), zlim=zlim.sd, xaxt="n", yaxt="n")
title(xlab="R-INLA", line=1)
image(sd2.spherical, col=viridis::inferno(100), zlim=zlim.sd, xaxt="n", yaxt="n")
title(xlab="Ovelapping RBMC", line=1)
dev.off()


t.idx = (1:n.t-1)*n.s+1
pdf("fig.sim.2.pdf", width=8, height=4)
par(mfrow=c(1,1), mar=c(2,2,1,1))
plot(abs(d2/d0-1)[t.idx], t="l", xlab=NA, ylab=NA, xaxt="n", yaxt="n")
axis(1, at=seq(0,n.t,length=5), labels=seq(0,n.t,length=5)*c(1,1,NA,1,1))
axis(2, at=seq(0,9e-4,length=2), labels=formatC(seq(0,9e-4,length=2), digits=0, format="e"))
title(xlab="Time", ylab="Relative error", line=1)
dev.off()

