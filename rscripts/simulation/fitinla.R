# command line arguments
args = commandArgs(trailingOnly=TRUE)
for(i in seq_along(args)){
    # problem dimensions
    if(args[i] == "-ms") m.s = as.integer(args[i+1])
    if(args[i] == "-mt") m.t = as.integer(args[i+1])
    
    # hyperparameters
    if(args[i] == "-rt") range.t = as.double(args[i+1])
    if(args[i] == "-rs") range.s = as.double(args[i+1])
    if(args[i] == "-sigma") sigma.st = as.double(args[i+1])
    if(args[i] == "-tauy") tau.y = as.double(args[i+1])
    if(args[i] == "-taub") tau.b = as.double(args[i+1])
    
    # number of cores
    if(args[i] == "-ncores1") n.cores1 = as.integer(args[i+1])
    if(args[i] == "-ncores2") n.cores2 = as.integer(args[i+1])
}

# load libraries
library(INLA)
library(INLAspacetime)
library(inlabru)
library(parallel)
inla.setOption(smtp = "pardiso", inla.mode = "compact", pardiso.license = "~/pardiso.license")
if(!exists(deparse(substitute(n.cores1)))) n.cores1 = detectCores()
if(!exists(deparse(substitute(n.cores2)))) n.cores2 = 1
source("generate.R")

# run INLA
data = data.frame(xcoord = rep(mesh.s$loc[,1], n.t),
                  ycoord = rep(mesh.s$loc[,2], n.t),
                  time = rep(1:n.t, each = n.s),
                  y = y, 
                  z = A.b[,2])
model = ~ -1 + Intercept(1) + z + field(list(space = cbind(xcoord, ycoord), time = time), model = stmodel)
stmodel = stModel.define(mesh.s, mesh.t, "121",
                         control.priors = list(prs = c(range.s, 0),
                                               prt = c(range.t, 0),
                                               psigma = c(sigma.st, 0)))
lkprec = list(prec = list(initial = log(tau.y), fixed = TRUE))
result = bru(model, 
             like(formula = y ~ ., 
                  family = "gaussian",
                  control.family = list(hyper = lkprec), 
                  data = data),
             options = list(verbose = TRUE,
                            safe = FALSE,
                            control.inla = list(int.strategy = "eb"),
                            control.fixed = list(prec = list(prec = 1e-5, prec.intercept = 1e-5))))

print(result$summary.fixed)
