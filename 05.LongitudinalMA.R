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
# LONGITUDINAL META-ANALYSIS (SECTION 4.5, TABLE 5 AND FIGURE 3)
################################################################################

# LOAD THE PACKAGE
library(mixmeta)

data(gliomas)
# THE gliomas DATASET IS ARRANGED IN A LONG FORMAT
head(gliomas)

# INDEPENDENT RANDOM EFFECTS, NO WITHIN-STUDY CORRELATION (MODEL 1)
mod1 <- mixmeta(logOR~0+factor(time), S=logORvar, random=~0+factor(time)|study,
  bscov="diag", data=gliomas)
print(summary(mod1), digits=3, report="var")

# COMPOUND-SYMMETRY RANDOM EFFECTS, NO WITHIN-STUDY CORRELATION (MODEL 2)
# NB: THIS REQUIRES A TWO-LEVEL MODEL WITH THE INNER-LEVEL VARIANCE FIXED TO 0
unit <- factor(seq(nrow(gliomas)))
mod2 <- mixmeta(logOR~0+factor(time), S=logORvar, random=~1|study/unit,
  bscov=c("unstr","fixed"), data=gliomas, control=list(Psifix=list(unit=0)))
print(summary(mod2), digits=3, report="var")

# HETEROGENEOUS AR1 RANDOM EFFECTS, NO WITHIN-STUDY CORRELATION (MODEL 3)
mod3 <- update(mod1, bscov="ar1")
print(summary(mod3), digits=3, report="var")

# BUILD THE HETEROGENEOUS AR1 WITHIN-STUDY ERRORS (CORRELATION AT 0.61)
cormat <- 0.61^abs(col(matrix(1,4,4)) - row(col(matrix(1,4,4))))
addS <- lapply(split(sqrt(gliomas$logORvar), gliomas$study), inputcov, cormat)
addS <- lapply(addS, function(x) x[apply(!is.na(x),1,any), 
  apply(!is.na(x),2,any)])

# INDEPENDENT RANDOM EFFECTS, HAR1 WITHIN-STUDY CORRELATION (MODEL 4)
mod4 <- mixmeta(logOR~0+factor(time), random=~0+factor(time)|study,
  bscov="diag", data=gliomas, control=list(addSlist=addS))
print(summary(mod4), digits=3, report="var")

# HAR1 RANDOM EFFECTS, HAR1 WITHIN-STUDY CORRELATION (MODEL 5)
mod5 <- update(mod4, bscov="ar1")
print(summary(mod5), digits=3, report="var")

# UNSTRUCTURED RANDOM EFFECTS, HAR1 WITHIN-STUDY CORRELATION (MODEL 6)
mod6 <- update(mod4, bscov="unstr")
print(summary(mod6), digits=3, report="var")

# COMPARE THE FIT WITH AIC
AIC(mod1, mod2, mod3, mod4, mod5, mod6)

# RE-RUN BEST FITTING MODEL WITH ML (ALLOWS TESTING OF FIXED EFFECTS)
mod4ml <- update(mod4, method="ml")
print(summary(mod4ml), digits=3, report="var")

# RANDOM-SLOPE MODEL WITH TIME AS CONTINUOUS AND CENTERED
mod7ml <- mixmeta(logOR~time, random=~I(time-15)|study, bscov="diag", 
  method="ml", data=gliomas, control=list(addSlist=addS, maxiter=200))
print(summary(mod7ml), digits=3, report="var")

# PREDICT
times <- unique(gliomas$time)
predmod4ml <- exp(predict(mod4ml, data.frame(time=times), ci=TRUE))
predmod7ml <- exp(predict(mod7ml, data.frame(time=times), ci=TRUE))

# PLOT
par(mar=c(5,4,1,0.5))
plot(c(0.5,2.5)~c(4,26), gliomas, type="n", xlab="Time (months)", 
  ylab="Survival OR")
abline(h=1)
colnew <- do.call(rgb,c(as.list(col2rgb("red")/255),list(0.2)))
polygon(c(times,rev(times)), c(predmod7ml[,2], rev(predmod7ml[,3])), col=colnew,
  border=NA)
arrows(times, predmod4ml[,2], times, predmod4ml[,3], col=4, angle=90, code=3, 
  length=0.05)
points(predmod4ml[,1]~times, pch=19, col=4)
lines(predmod7ml[,1]~times, type="o", pch=19, col=2)
legend("top",c("Model 4 (Indicators)","Model 7 (continuous)"), col=c(4,2),
  lty=c(NA,1),pch=19, cex=0.8, ncol=2, bty = "n", inset=0.05)

# WE COMPARE THE TWO MODELS
AIC(mod4ml, mod7ml)

# SEE help(gliomas) AND help(dbs) FOR FURTHER INFO
