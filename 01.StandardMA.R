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
# STANDARD META-ANALYSIS (SECTION 4.1)
################################################################################

# LOAD THE PACKAGES
library(mixmeta); library(Epi)

# STANDARD RANDOM-EFFECTS META-ANALYSIS WITH MAXIMUM LIKELIHOOD
uniran <- mixmeta(logor, logorvar, data=bcg, method="ml")

# RESULTS
print(summary(uniran), digits=3, report="var")
print(ci.exp(uniran), digits=3)

# EXTRACT LOG-OR AND CALCULATE 95% CONFIDENCE INTERVALS
pred <- with(bcg, cbind(logor, logor-1.96*sqrt(logorvar),
  logor+1.96*sqrt(logorvar)))

# BEST-LINEAR UNBIASED PREDICTIONS, WITH PREDICTION INTERVALS 
blup <- blup(uniran, pi=TRUE)

# FOREST PLOT
plot(pred[,1], rev(bcg$trial)+0.2, xlim=c(-3,3), ylim=c(0,14), pch=18,
  axes=FALSE, xlab="Log odds ratio", ylab="Trial", main="Forest plot")
axis(1)
axis(2, at=bcg$trial, labels=rev(bcg$trial), lty=0, las=1)
abline(v=coef(uniran))
segments(pred[,2], rev(bcg$trial)+0.2, pred[,3], rev(bcg$trial)+0.2, lty=5)
points(blup[,1], rev(bcg$trial)-0.2, pch=19)
segments(blup[,2], rev(bcg$trial)-0.2, blup[,3], rev(bcg$trial)-0.2)
legend("right", c("Original","BLUPs"),lty=c(5,1), pch=c(18,19), lwd=1, bty="n")

# META-REGRESSION MODEL TO EVALUATE THE EFFECT OF LATITUDE
uniranlat <- update(uniran, .~. + ablat)

# LIKELIHOOD RATIO TEST (ALLOWED WITH ML)
drop1(uniranlat, test="Chisq")

# RESULTS
print(summary(uniranlat), digits=3, report="var")

# SEE help(bcg) FOR FURTHER INFO
