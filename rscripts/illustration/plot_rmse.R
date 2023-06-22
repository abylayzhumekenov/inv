# suppressWarnings(suppressMessages(library(INLA)))

n = 99
phi = 0.95

A = 1:49
S = 50
B = 51:99
AA = 1:59
SA = 60
BB = 41:99
SB = 40

samples = 10^(1:3)
cols = hcl.colors(3, rev=TRUE)
cols_ = hcl.colors(3, rev=TRUE, alpha=0.5)

pdf("img/fig.var.4.pdf", width = 14, height = 5)
par(mfrow=c(1,2), mar=c(3,3,1,1))

# plot 1
rmse1 = phi^(2*abs((1:n)-S))
plot(rmse1, xlim=c(1,n), ylim=c(0,1), t="l", lwd=2, xlab=NA, ylab=NA)
for(k in seq_along(samples)){
    rmse2 = rep(sqrt(2/samples[k]), length(rmse1))
    lines(rmse2, lwd=2, lty=3, col=cols_[k])
}
for(k in seq_along(samples)){
    rmse3 = rmse1*sqrt(2/samples[k])
    lines(rmse3, lwd=2, lty=1, col=cols[k])
}
title(xlab="Index", ylab="Relative RMSE", line=2)
legend("topright", 
       legend=c("MC", "RBMC, no samples", "RBMC, 10 samples", "RBMC, 100 samples", "RBMC, 1000 samples"), 
       col=c("gray", "black", hcl.colors(3, rev=TRUE)), bty="n",
       lwd=c(2,2,2,2,2), lty=c(3,1,1,1,1))

# plot 2
rmse1_AA = phi^(2*abs(1:S-SA))
rmse1_BB = phi^(2*abs(S:n-SB))
rmse1_SA = phi^(2*abs(S:SA-SA))
rmse1_SB = phi^(2*abs(SB:S-SB))
plot(1:S, rmse1_AA, xlim=c(1,n), ylim=c(0,1), t="l", lwd=2, xlab=NA, ylab=NA)
lines(S:n, rmse1_BB, lwd=2, col=rgb(0,0,0,1))
lines(S:SA, rmse1_SA, lwd=2, lty=2, col=rgb(0,0,0,0.5))
lines(SB:S, rmse1_SB, lwd=2, lty=2, col=rgb(0,0,0,0.5))
for(k in seq_along(samples)){
    rmse2 = rep(sqrt(2/samples[k]), length(rmse1))
    lines(rmse2, lwd=2, lty=3, col=cols_[k])
}
for(k in seq_along(samples)){
    rmse3_AA = rmse1_AA*sqrt(2/samples[k])
    rmse3_BB = rmse1_BB*sqrt(2/samples[k])
    rmse3_SA = rmse1_SA*sqrt(2/samples[k])
    rmse3_SB = rmse1_SB*sqrt(2/samples[k])
    lines(1:S, rmse3_AA, lwd=2, lty=1, col=cols[k])
    lines(S:n, rmse3_BB, lwd=2, lty=1, col=cols[k])
    lines(S:SA, rmse3_SA, lwd=2, lty=2, col=cols_[k])
    lines(SB:S, rmse3_SB, lwd=2, lty=2, col=cols_[k])
}
title(xlab="Index", ylab="Relative RMSE", line=2)
legend("topright", 
       legend=c("MC", "Discarded", "RBMC, no samples", "RBMC, 10 samples", "RBMC, 100 samples", "RBMC, 1000 samples"), 
       col=c("gray", rgb(0,0,0,0.5), "black", hcl.colors(3, rev=TRUE)), bty="n",
       lwd=c(2,2,2,2,2,2), lty=c(3,2,1,1,1,1))
dev.off()
