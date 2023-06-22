pdf("img/fig.var.2.pdf", width=10, height=5)
par(mfrow=c(1,2), mar=c(1,1,1,1))

plot(0, col=NA, xlim=c(0,11), ylim=c(0,11), asp=1, xlab=NA, ylab=NA, axes=FALSE)
polygon(c(5,6,6,5), c(0,0,11,11), col=rgb(0,0,0,1), border=NA)
polygon(c(0,11,11,0), c(0,0,11,11), lwd=2)
text(5.5, 5.5, "S", col="white")
text(c(2.5,8.5), c(5.5,5.5), c(expression(A[1]),expression(A[2])))
# text(5.5, 5.5, expression(x[S]), col="white")
# text(c(2.5,8.5), c(5.5,5.5), c(expression(x[A]), expression(x[B])))

plot(0, col=NA, xlim=c(0,11), ylim=c(0,11), asp=1, xlab=NA, ylab=NA, axes=FALSE)
polygon(c(0,11,11,0), c(5,5,6,6), col=rgb(0,0,0,1), border=NA)
polygon(c(5,6,6,5), c(0,0,11,11), col=rgb(0,0,0,1), border=NA)
polygon(c(0,11,11,0), c(0,0,11,11), lwd=2)
text(5.5, 5.5, "S", col="white")
text(c(2.5,2.5,8.5,8.5), c(8.5,2.5,8.5,2.5), c(expression(A[1]),expression(A[2]),expression(A[3]),expression(A[4])))

dev.off()