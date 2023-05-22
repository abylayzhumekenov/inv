# ------------------------------------------------------------------------------
# ---------- LOAD THE RESULTS  -------------------------------------------------
# ------------------------------------------------------------------------------

# load the data
suppressMessages(suppressWarnings(library(INLA)))
load("data/data.Rdata")
load("data/theta.Rdata")
load("data/result.Rdata")
if(TRUE){           # rewrite?
    mu.inv = readBin("../../data/mu", "double", n.s*n.t+n.b+1, endian = "swap")[-1]
    d.inv = readBin("../../data/d", "double", n.s*n.t+n.b+1, endian = "swap")[-1]
    sd.inv = sqrt(d.inv)
    save(list = c("mu.inla", "sd.inla", "mu.inv", "sd.inv"), file = "data/result.Rdata")
}

# ------------------------------------------------------------------------------
# ---------- CREATE PROJECTIONS  -----------------------------------------------
# ------------------------------------------------------------------------------

# create a projection
plot.res = c(256, 128)
plot.asp = plot.res[2] / plot.res[1]
coords.spherical = expand.grid(seq(0,2*pi,length=plot.res[1]+1)[-1], seq(0,pi,length=plot.res[2]+1)[-1])
coords.cartesian = data.frame(x = sin(coords.spherical$Var2) * cos(coords.spherical$Var1),
                              y = sin(coords.spherical$Var2) * sin(coords.spherical$Var1),
                              z = cos(coords.spherical$Var2))
A.spherical = inla.spde.make.A(mesh.s, as.matrix(coords.cartesian))

# project the mean
mu.inla.spherical = drop(A.spherical %*% mu.inla[1:n.s])
mu.inla.spherical = matrix(mu.inla.spherical, plot.res[1])
mu.inv.spherical = drop(A.spherical %*% mu.inv[1:n.s])
mu.inv.spherical = matrix(mu.inv.spherical, plot.res[1])

# project the sd
sd.inla.spherical = drop(A.spherical %*% sd.inla[1:n.s])
sd.inla.spherical = matrix(sd.inla.spherical, plot.res[1])
sd.inv.spherical = drop(A.spherical %*% sd.inv[1:n.s])
sd.inv.spherical = matrix(sd.inv.spherical, plot.res[1])

# set the plot range
zlim.mu = range(mu.inla.spherical, mu.inv.spherical)
zlim.sd = range(sd.inla.spherical, sd.inv.spherical)

# ------------------------------------------------------------------------------
# ---------- PLOT THE RESULTS  -------------------------------------------------
# ------------------------------------------------------------------------------

# plot the spatial field
pdf("img/fig.sim.1.pdf", width=9, height=5)
par(mfrow=c(2,2), mar=c(2,2,1,1))
image(mu.inla.spherical, col=viridis::viridis(100), asp=plot.asp, zlim=zlim.mu, xaxt="n", yaxt="n", axes=FALSE)
title(ylab="Mean", line=1)
image(mu.inv.spherical, col=viridis::viridis(100), asp=plot.asp, zlim=zlim.mu, xaxt="n", yaxt="n", axes=FALSE)
image(sd.inla.spherical, col=viridis::inferno(100), asp=plot.asp, zlim=zlim.sd, xaxt="n", yaxt="n", axes=FALSE)
title(xlab="R-INLA", ylab="SD", line=1)
image(sd.inv.spherical, col=viridis::inferno(100), asp=plot.asp, zlim=zlim.sd, xaxt="n", yaxt="n", axes=FALSE)
title(xlab="Ovelapping RBMC", line=1)
dev.off()

# plot the temporal field
pdf("img/fig.sim.2.pdf", width=5, height=5)
par(mfrow=c(1,1), mar=c(2,2,1,1))
t.idx = (1:n.t-1)*n.s+1
error = abs(sd.inv/sd.inla-1)[t.idx]
plot(error, t="l", col=rgb(0.2,0.4,1,1), lwd=2, xlab=NA, ylab=NA, xaxt="n", yaxt="n", axes=FALSE)
axis(1, at=seq(0,n.t,length=5), labels=seq(0,n.t,length=5)*c(1,1,NA,1,1), col="gray")
axis(2, at=seq(0,max(error),length=2), labels=c(0, formatC(max(error), digits=0, format="e")), col="gray")
title(xlab="Time", ylab="Relative error", line=1)
dev.off()



# plot the strong and weak scaling
pdf("img/fig.sim.3.pdf", width=11, height=5)
par(mfrow=c(1,2), mar=c(3,3,1,1))

p1 = 2^(1:10)
t1 = c(588.62, 331.84, 194.09, 89.06, 51.12, 29.53, 19.73, 13.40, 11.16, 13.30)
plot(log2(p1), log2(p1/p1[1]), t="l", lwd=2, lty=3, col=rgb(1,0.2,0.2,1), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
lines(log2(p1), log2(t1[1]/t1), lwd=2, col=rgb(0.2,0.4,1,1))
points(log2(p1), log2(t1[1]/t1), pch=16, cex=1, col=rgb(0.2,0.4,1,1))
axis(1, at=log2(p1), labels=p1, line=-0.5, tick=FALSE)
axis(2, at=log2(p1/p1[1]), labels=p1/p1[1], line=-0.5, tick=FALSE)
title(xlab="Processors", ylab="Speedup", line=2)
legend(x=5, y=3, legend=c("Ideal","Real"), pch=c(NA,16), lty=c(3,1), lwd=c(2,2), col=c(rgb(1,0.2,0.2,1),rgb(0.2,0.4,1,1)), bty="n")

p2 = 2^(1:10)
t2 = c(11.71, 12.42, 13.50, 15.53, 19.74, 20.83, 24.56, 31.97, 35.34, 49.10)
plot(log2(p2), rep(1, length(p2)), ylim=c(0,1), t="l", lwd=2, lty=3, col=rgb(1,0.2,0.2,1), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
lines(log2(p2), t2[1]/t2, lwd=2, col=rgb(0.2,0.4,1,1))
points(log2(p2), t2[1]/t2, pch=16, cex=1, col=rgb(0.2,0.4,1,1))
axis(1, at=log2(p2), labels=p2, line=-0.5, tick=FALSE)
axis(2, at=seq(0,1,0.25), labels=seq(0,1,0.25), line=-0.5, tick=FALSE)
title(xlab="Processors", ylab="Efficiency", line=2)
legend(x=2, y=0.2, legend=c("Ideal","Real"), pch=c(NA,16), lty=c(3,1), lwd=c(2,2), col=c(rgb(1,0.2,0.2,1),rgb(0.2,0.4,1,1)), bty="n")

dev.off()
