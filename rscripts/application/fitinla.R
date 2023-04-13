# command line arguments
args = commandArgs(trailingOnly=TRUE)
for(i in seq_along(args)){
    # problem dimensions
    if(args[i] == "-ms") m.s = as.integer(args[i+1])
    if(args[i] == "-mt") m.t = as.integer(args[i+1])
    if(args[i] == "-res") res.s = as.integer(args[i+1])
    if(args[i] == "-nh") n.h = as.integer(args[i+1])
    
    # hyperparameters
    if(args[i] == "-rt") range.t = as.double(args[i+1])
    if(args[i] == "-rs") range.s = as.double(args[i+1])
    if(args[i] == "-sigma") sigma = as.double(args[i+1])
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

# prepare the data
source("prepdata.R")

# define a model
model = ~ -1 + Intercept(1) + elevation + northing
for(i in 1:n.h) model = update(model, paste("~ . +", paste0("harmonic", i, ".sin"), " + ", paste0("harmonic", i, ".cos")))
model = update(model, ~ . + field(list(space = cbind(xcoord, ycoord), time = time), model = stmodel))

# set hyperparameters
if(!exists(deparse(substitute(r.s)))) r.s = 1259
if(!exists(deparse(substitute(r.t)))) r.t = 50
if(!exists(deparse(substitute(sigma.st)))) sigma.st = 5.96
if(!exists(deparse(substitute(tau.y)))) tau.y = 0.143

# define priors
stmodel = stModel.define(mesh.s, mesh.t, "121",
                         control.priors = list(prs = c(r.s, 0),
                                               prt = c(r.t, 0),
                                               psigma = c(sigma.st, 0)))
lkprec = list(prec = list(initial = log(tau.y), fixed = TRUE))

# fit the model
result = bru(model, 
             like(formula = y ~ ., 
                  family = "gaussian",
                  control.family = list(hyper = lkprec), 
                  data = data),
             options = list(verbose = TRUE,
                            safe = FALSE,
                            num.threads = paste0(n.cores1, ":", n.cores2),
                            control.inla = list(cmin=0, int.strategy = "eb")))

# print the result
print(result$summary.fixed)
print(result$internal.summary.hyperpar)
print(result$cpu.used)
