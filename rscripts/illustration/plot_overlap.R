pdf("img/fig.var.2.pdf", width=10, height=5)
par(mfrow=c(1,2), mar=c(1,1,1,1))

# plot(0, col=NA, xlim=c(0,11), ylim=c(0,11), asp=1, xlab=NA, ylab=NA, axes=FALSE)
# polygon(c(5,6,6,5), c(0,0,11,11), col=rgb(0,0,0,1), border=NA)
# polygon(c(0,11,11,0), c(5,5,6,6), col=rgb(0,0,0,1), border=NA)
# polygon(c(0,11,11,0), c(0,0,11,11), lwd=2)
# text(5.5, 5.5, "S", col="white")
# text(c(2.5,2.5,8.5,8.5), c(8.5,2.5,8.5,2.5), c(expression(A[1]), expression(A[2]), expression(A[3]), expression(A[4])))

plot(0, col=NA, xlim=c(0,11), ylim=c(0,11), asp=1, xlab=NA, ylab=NA, axes=FALSE)
polygon(c(5,6,6,5), c(0,0,4,4), col=rgb(0.5,0.5,0.5,1), border=NA)
polygon(c(7,11,11,7), c(5,5,6,6), col=rgb(0.5,0.5,0.5,1), border=NA)
polygon(c(0,7,7,6,6,0), c(4,4,11,11,5,5), col=rgb(0,0,0,1), border=NA)
# polygon(c(0,5,5,0), c(6,6,11,11), col=NA, border=rgb(0.5,0.5,0.5,1), lwd=1.5, lty=3)
polygon(c(0,11,11,0), c(0,0,11,11), lwd=2)
text(6.5, 4.5, expression(S[1]), col="white")
text(c(3,2.5,8.5,8.5), c(8,2.5,8.5,2.5), c(expression(A[1]), expression(A[2]), expression(A[3]), expression(A[4])),
     col=c(rgb(0,0,0,1), rgb(0.5,0.5,0.5,1), rgb(0.5,0.5,0.5,1), rgb(0.5,0.5,0.5,1)))

plot(0, col=NA, xlim=c(0,11), ylim=c(0,11), asp=1, xlab=NA, ylab=NA, axes=FALSE)
polygon(c(5,6,6,5), c(7,7,11,11), col=rgb(0.5,0.5,0.5,1), border=NA)
polygon(c(7,11,11,7), c(5,5,6,6), col=rgb(0.5,0.5,0.5,1), border=NA)
polygon(c(0,6,6,7,7,0), c(6,6,0,0,7,7), col=rgb(0,0,0,1), border=NA)
polygon(c(0,11,11,0), c(0,0,11,11), lwd=2)
text(6.5, 6.5, expression(S[2]), col="white")
text(c(2.5,3,8.5,8.5), c(8.5,3,8.5,2.5), c(expression(A[1]), expression(A[2]), expression(A[3]), expression(A[4])),
     col=c(rgb(0.5,0.5,0.5,1), rgb(0,0,0,1), rgb(0.5,0.5,0.5,1), rgb(0.5,0.5,0.5,1)))

dev.off()
