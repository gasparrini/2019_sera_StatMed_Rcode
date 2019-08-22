################################################################################
# Updated R code for the analysis in:
#
#   "An extended mixed-effects framework for meta-analysis"
#   Francesco Sera, Ben Armstrong & Antonio Gasparrini
#   Statistics in Medicine - 2019
#   http://www.ag-myresearch.com/2019_sera_statmed.html
#
# Update: 22 Aug 2019
# * an up-to-date version of this code is available at:
#   https://github.com/gasparrini/2019_sera_StatMed_Rcode
################################################################################

################################################################################
# SIMULATION STUDY (SECTION 6, TABLE 6)
# WITH PARALLELIZATION: THIS ROUTINE MAY NOT WORK IN SOME OS AND MACHINE
################################################################################

# LOAD THE PACKAGES
library(mixmeta) ; library(foreach) ; library(doParallel)
library(iterators) ; library(parallel)

# DEFINE FIXED PARAMETERS
g2 <- 10
beta <- 0
tau <- 1

# DATA FRAME WITH COMBINATIONS OF SCENARIOS
comb <- data.frame(m=c(10,50), rhob=rep(c(0,0.8),each=2),
  rhow=rep(c(0,0.8),each=4))
comb$rhob <- ifelse(comb$rhob <= -1/(3-1), -1/3, comb$rhob)
comb$rhow <- ifelse(comb$rhow <= -1/(3-1), -1/3, comb$rhow)

# NAMES
stats <- c("bias","rmse","cov")
ran <- c("tau1","tau2","rhob1","rhob2")

# NOMINAL VALUE
qn <- qnorm(0.975)

# NUMBER OF SIMULATIONS
nsim <- 10000

################################################################################

# PREPARE THE PARALLELIZATION
ncores <- detectCores()
cl <- makeCluster(max(1,ncores-2))
registerDoParallel(cl)

################################################################################

# START THE NESTED LOOP BY SCENARIO/ITERATIONS
temp <- foreach(combi=iter(comb,by="row"), .packages=c("mixmeta")) %:% 
  foreach(i=icount(nsim), .combine=rbind) %dopar% {

  # DEFINE THE DATA FROM WHICH TO SIMULATE
  n <- combi$m * g2
  y <- matrix(0,n,3)
  S <- inputcov(matrix(runif(n*3, 0.5, 2), n, 3), cor=combi$rhow)
  Psi1 <- Psi2 <- inputcov(rep(tau, 3), cor=combi$rhob)
  level1 <- rep(seq(combi$m), each=g2)
  level2 <- rep(seq(g2), combi$m)
  
  # SIMULATE (WITH SEED)
  set.seed(13041975+i)
  sim <- mixmetaSim(y,S,Psi=list(Psi1,Psi2),random=~1|level1/level2)
  
  # FIT THE MODEL (PREVENT ERRORS DUE TO NON-CONVERGENCE)
  arglist <- list(sim~1, S=S, random=~1|level1/level2, bscov="cs")
  model <- tryCatch(do.call("mixmeta",arglist), error=function(x) NULL)
  
  # RETURN THE VECTOR OF ESTIMATES (SET TO NA IF NON-CONVERGENCE)
  if(!is.null(model)) c(coef(model)[1], sqrt(vcov(model)[1,1]), 
    model$Psi[[1]][1,1], model$Psi[[2]][1,1], cov2cor(model$Psi[[1]])[1,2],
    cov2cor(model$Psi[[2]])[1,2]) else NA
}

################################################################################

# REMOVE PARALLELIZATION
stopCluster(cl)

################################################################################

# COMPUTE THE STATS (REMOVING THE MISSINGS)
res <- t(sapply(seq(temp), function(i) {
  nan <- sum(!is.na(temp[[i]][,1]))
  x <- na.omit(temp[[i]])
  bias <- colMeans(x[,-2]) - c(beta,tau,tau,rep(comb[i,"rhob"],2))
  c(bias[1],
    sqrt(mean((x[,1]-beta)^2)), 
    mean(x[,1]-qn*x[,2]<=beta & x[,1]+qn*x[,2]>=beta),
    bias[-1],
    nan/nsim)
}))
colnames(res) <- c(paste("beta",stats,sep="-"),paste(ran,"bias",sep="-"),"conv")

################################################################################

# SAVE THE RESULTS
#save.image("simul.RData")
