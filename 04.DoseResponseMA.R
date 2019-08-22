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
# DOSE-RESPONSE META-ANALYSIS (SECTION 4.4, TABLE 4 AND FIGURE 2)
################################################################################

# LOAD THE PACKAGES
library(mixmeta); library(dosresmeta); library(splines)

# INSPECT THE DATA
head(alcohol)

# COMPUTE THE WITHIN-STUDY CORRELATIONS EXCLUDING THE REFERENCE
addS <- lapply(split(alcohol, alcohol$id), function(x)
  covar.logrr(y=logrr, v=se^2, cases=cases, n=peryears, type=type, data=x))
sub <- subset(alcohol, !is.na(se))

# LINEAR FIXED AND RANDOM EFFECTS NOT ACCOUNTING FOR WITHIN-STUDY CORRELATIONS 
modL1 <- mixmeta(logrr ~ 0 + dose, S=se^2, random= ~ 0 + dose|id, data=sub,
  method="ml")
summary(modL1)

# LINEAR FIXED AND RANDOM EFFECTS ACCOUNTING FOR WITHIN-STUDY CORRELATIONS
modL2 <- mixmeta(logrr ~ 0 + dose, random= ~ 0 + dose|id, data=sub, method="ml",
  control=list(addSlist=addS))
summary(modL2)

# NON-LINEAR FIXED AND RANDOM EFFECTS
modNL1 <- mixmeta(logrr ~ 0 + ns(dose, knots=c(10,25)), data=sub, 
  random= ~ 0 + ns(dose, knots=c(10,25))|id, method="ml",
  control=list(addSlist=addS))
summary(modNL1)

# SIMPLIFY THE MODEL BY ALLOWING NON-LINEARITY ONLY IN FIXED EFFECTS
modNL2 <- update(modNL1, random= ~ 0 + dose|id)
summary(modNL2)

# FIXED-EFFECTS MODEL (TRICK: random TO DEFINE THE GROUPING, THEN FIX IT TO 0)
modNL3 <- mixmeta(logrr ~ 0 + ns(dose, knots=c(10,25)), random= ~ 1|id,
  data=sub, method="ml",bscov="fixed", control=list(addSlist=addS, Psifix=0))
summary(modNL3)

# COMPARE WITH AIC
AIC(modL1, modL2, modNL1, modNL2, modNL3)

# PREDICT THE RR FOR 12g/day FOM TWO MODELS
exp(predict(modL1, newdata=data.frame(dose=12), ci=TRUE))
exp(predict(modL2, newdata=data.frame(dose=12), ci=TRUE))

# PREDICT THE RR ALONG THE DOSE RANGE
predlin <- exp(predict(modL2, newdata=data.frame(dose=0:60), ci=TRUE))
prednonlin <- exp(predict(modNL1, newdata=data.frame(dose=0:60), ci=TRUE))

# PLOT
par(mar=c(5,4,1,0.5))
col1 <- do.call(rgb,c(as.list(col2rgb("blue")/255), list(0.2)))
col2 <- do.call(rgb,c(as.list(col2rgb("green")/255), list(0.2)))
plot(0:60,predlin[,1], type="l", ylim=c(0.85,1.9), ylab="RR",
  xlab="Alcohol intake (gr/day)")
polygon(c(0:60,60:0), c(predlin[,2],rev(predlin[,3])), col=col1 ,border=NA)
lines(0:60,prednonlin[,1], lty=5)
polygon(c(0:60,60:0), c(prednonlin[,2], rev(prednonlin[,3])), col=col2 ,border=NA)
legend("topleft", c("Model L2","Model NL1"), lty= c(1,5), bty="n", inset=0.1)

# SEE help(alcohol) FOR FURTHER INFO
