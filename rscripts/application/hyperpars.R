library(geoR)
source("generate.R")

y[id.na] = 0
w = drop(solve(crossprod(A.b, A.b), crossprod(A.b, y)))
yhat = drop(A.b %*% w)
r = y - yhat

nvg.t = 20
nvg.s = 1000
id.t = sample(1:m.t, nvg.t)
id.s = sample(1:m.s, nvg.s)
id.st = rep(id.t * m.s, each = nvg.s) + rep(id.s, nvg.t)
loc.vg = cbind(rep(loc[id.s,1], nvg.t), rep(loc[id.s,2], nvg.t))
colnames(loc.vg) = colnames(loc)
r.vg = r[id.st]
vg = variog(coords = loc.vg, data = r.vg)

# plot(vg)