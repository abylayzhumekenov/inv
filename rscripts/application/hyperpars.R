library(geoR)
source("generate.R")

AA.b = cbind(x0, x1, x2, x3)
AA.b[id.na,] = 0
yy = y
yy[id.na] = 0
w = drop(solve(crossprod(AA.b, AA.b), crossprod(AA.b, yy)))
yhat = drop(A.b %*% w)
r = y - yhat

# sigma.sq
arfit = ar(r[(1:m.t-1)*m.s+1000], order.max = 1, na.action = na.pass)
tau_y = 1/arfit$var.pred
sigma.sq = var(r[(1:m.t-1)*m.s+1000], na.rm = TRUE) - arfit$var.pred

# range.t
# plot(acf(r[(1:m.t-1)*m.s+10], m.t, na.action = na.pass))
range.t = 2

# range.s
nvg.t = 20
nvg.s = 1000
id.t = sample(1:m.t, nvg.t)
id.s = sample(1:m.s, nvg.s)
id.st = rep(id.t * m.s, each = nvg.s) + rep(id.s, nvg.t)
loc.vg = cbind(rep(loc[id.s,1], nvg.t), rep(loc[id.s,2], nvg.t))
colnames(loc.vg) = colnames(loc)
r.vg = r[id.st]
loc.vg = loc.vg[!is.na(r.vg),]
r.vg = r.vg[!is.na(r.vg)]
vg = variog(coords = loc.vg, data = r.vg)
# plot(vg)
range.s = 1500.000

print(c(range.t, range.s, sigma.sq, tau_y, tau_b))
source("generate.R")
