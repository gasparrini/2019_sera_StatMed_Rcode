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
################################################################################

# LOAD THE PACKAGES
library(mixmeta)

# DEFINE FIXED PARAMETERS
g2 <- 10
beta <- 0
tau <- 1

# DATA FRAME WITH COMBINATIONS OF SCENARIOS
comb <- data.frame(m=c(10,50), rhob=rep(c(0,0.8),each=2),
  rhow=rep(c(0,0.8),each=4))
comb$rhob <- ifelse(comb$rhob <= -1/(3-1), -1/3, comb$rhob)
comb$rhow <- ifelse(comb$rhow <= -1/(3-1), -1/3, comb$rhow)

# MATRIX WITH FINAL RESULTS
res <- matrix(NA, nrow(comb), 8)
stats <- c("bias","rmse","cov")
ran <- c("tau1","tau2","rhob1","rhob2")
colnames(res) <- c(paste("beta",stats,sep="-"),paste(ran,"bias",sep="-"),"conv")

# NOMINAL VALUE
qn <- qnorm(0.975)

# NUMBER OF SIMULATIONS
nsim <- 1000

################################################################################

# START THE LOOP BY SCENARIO
for(i in seq(nrow(comb))) {
  
  # PRINT
  cat("\n\n ",paste("Combination",i),"\n")
  
  # DEFINE THE DATA FROM WHICH TO SIMULATE
  n <- comb$m[i] * g2
  y <- matrix(0,n,3)
  S <- inputcov(matrix(runif(n*3, 0.5, 2), n, 3), cor=comb$rhow[i])
  Psi1 <- Psi2 <- inputcov(rep(tau, 3), cor=comb$rhob[i])
  level1 <- rep(seq(comb$m[i]), each=g2)
  level2 <- rep(seq(g2), comb$m[i])
  
  # SIMULATE (WITH SEED)
  set.seed(13041975+i)
  simlist <- mixmetaSim(y,S,Psi=list(Psi1,Psi2),random=~1|level1/level2,nsim=nsim)
  
  # BUILD THE TEMP0RARY OBJECT TO STORE THE ESTIMATES FROM EACH MODEL
  temp <- matrix(NA, nsim, 6, dimnames=list(NULL,c("beta","beta.se",ran)))

################################################################################
  
  # START THE LOOP BY ITERATION
  for(j in seq(nsim)) {
    
    # PRINT
    cat(j,"")

    # FIT THE MODEL (PREVENT ERRORS DUE TO NON-CONVERGENCE)
    arglist <- list(simlist[[j]]~1, S=S, random=~1|level1/level2, bscov="cs")
    model <- tryCatch(do.call("mixmeta",arglist), error=function(x) NULL)
    
    # STORE THE ESTIMATES (SET TO NA IF NON-CONVERGENCE)
    temp[j,] <- if(!is.null(model)) c(coef(model)[1], sqrt(vcov(model)[1,1]), 
      model$Psi[[1]][1,1], model$Psi[[2]][1,1], cov2cor(model$Psi[[1]])[1,2],
      cov2cor(model$Psi[[2]])[1,2]) else NA
  }
  
################################################################################
  
  # COMPUTE THE STATS (REMOVING THE MISSINGS)
  res[i,8] <- sum(!is.na(temp[,1]))/nsim
  temp2 <- na.omit(temp)
  res[i,-c(2,3,8)] <- colMeans(temp2[,-2]) - c(beta,tau,tau,rep(comb$rhob[i],2))
  res[i,2] <- sqrt(mean((temp2[,1]-beta)^2))
  res[i,3] <- mean(temp2[,1]-qn*temp2[,2]<=beta & temp2[,1]+qn*temp2[,2]>=beta)
  
}

################################################################################

# SAVE THE RESULTS
#save.image("simul.RData")
