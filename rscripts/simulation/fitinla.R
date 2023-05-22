# ------------------------------------------------------------------------------
# ---------- SET PARAMETERS ----------------------------------------------------
# ------------------------------------------------------------------------------

# command line arguments
args = commandArgs(trailingOnly=TRUE)
for(i in seq_along(args)){
    # fixed hyperparameters
    if(args[i] == "-prs") prs = as.double(args[i+1])
    if(args[i] == "-prt") prt = as.double(args[i+1])
    if(args[i] == "-psig") psig = as.double(args[i+1])
    if(args[i] == "-ptauy") ptauy = as.logical(as.integer(args[i+1]))
    
    # INLA parameters
    if(args[i] == "-v") verbose = as.logical(as.integer(args[i+1]))
    if(args[i] == "-ncores1") n.cores1 = as.integer(args[i+1])
    if(args[i] == "-ncores2") n.cores2 = as.integer(args[i+1])
}

if(!exists(deparse(substitute(prs)))) prs = 0
if(!exists(deparse(substitute(prt)))) prt = 0
if(!exists(deparse(substitute(psig)))) psig = 0
if(!exists(deparse(substitute(ptauy)))) ptauy = FALSE

# ------------------------------------------------------------------------------
# ---------- LOAD LIBRARIES AND DATA -------------------------------------------
# ------------------------------------------------------------------------------

# load libraries
suppressMessages(suppressWarnings(library(INLA)))
library(INLAspacetime)
library(inlabru)
library(parallel)
inla.setOption(smtp = "pardiso", inla.mode = "compact", pardiso.license = "~/pardiso.license")
if(!exists(deparse(substitute(n.cores1)))) n.cores1 = detectCores()
if(!exists(deparse(substitute(n.cores2)))) n.cores2 = 1
if(!exists(deparse(substitute(verbose)))) verbose = FALSE

# load the data
load("data/data.Rdata")

# ------------------------------------------------------------------------------
# ---------- FIT THE MODEL USING INLA ------------------------------------------
# ------------------------------------------------------------------------------

# prepare the data
data = data.frame(xcoord = rep(mesh.s$loc[,1], n.t),
                  ycoord = rep(mesh.s$loc[,2], n.t),
                  zcoord = rep(mesh.s$loc[,3], n.t),
                  time = rep(1:n.t, each = n.s),
                  y = y)

# define the model
model = ~ -1 + Intercept(1) + field(list(space = cbind(xcoord, ycoord, zcoord), time = time), model = stmodel)
stmodel = stModel.define(mesh.s, mesh.t, "121",
                         control.priors = list(prs = c(r.s, prs),
                                               prt = c(r.t, prt),
                                               psigma = c(sig, psig)))
lkprec = list(prec = list(initial = log(tau.y), fixed = !ptauy))

# run inla
result = bru(model, 
             like(formula = y ~ ., 
                  family = "gaussian",
                  control.family = list(hyper = lkprec), 
                  data = data),
             options = list(verbose = verbose,
                            safe = FALSE,
                            num.threads = paste0(n.cores1, ":", n.cores2),
                            control.inla = list(int.strategy = "eb"),
                            control.fixed = list(prec = list(prec = 1e-5, prec.intercept = 1e-5))))

# ------------------------------------------------------------------------------
# ---------- SAVE THE RESULTS  -------------------------------------------------
# ------------------------------------------------------------------------------

# save the posterior latent field
mu.inla = c(result$summary.random$field$mode, result$summary.fixed$mode)
sd.inla = c(result$summary.random$field$sd, result$summary.fixed$sd)
save(list = c("mu.inla", "sd.inla"), file = "data/result.Rdata")

# save the posterior hyperparameters
theta.inla = c(r.s,r.t,sig,tau.y)
if(!is.null(result$internal.summary.hyperpar)){
    theta.inla[as.logical(c(prs,prt,psig,ptauy))] = exp(result$internal.summary.hyperpar$mode)
}
save(list = c("theta.inla"), file = "data/theta.Rdata")

# print out info
cat(format(c("ns","nt","nu","ms","mt","mu","nb"), width=12, justify="right"), "\n",
    format(c(n.s,n.t,n.s*n.t,m.s,m.t,m.s*m.t,n.b), width=12, justify="right"), "\n",
    "\n",
    format(c("rs","rt","sig","tauy","taub"), width=12, justify="right"), "\n",
    format(c(theta.inla,tau.b), width=12, justify="right"), "\n",
    "\n",
    format(names(result$cpu.used), width=12, justify="right"), "\n",
    format(result$cpu.used, width=12, justify="right"), "\n",
    sep = "")

# print(head(diag(result$misc$configs$config[[1]]$Qprior)))