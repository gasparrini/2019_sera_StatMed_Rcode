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
# MULTIVARIATE (NETWORK) META-ANALYSIS (SECTION 4.2, TABLE 1)
################################################################################

# LOAD THE PACKAGE
library(mixmeta)

# INSPECT THE DATA
head(smoking)
names(smoking)

# CONSISTENCY MODEL, UNSTRUCTURED BETWEEN-STUDY (CO)VARIANCE
y <- as.matrix(smoking[11:13])
S <- as.matrix(smoking[14:19])
mod1 <- mixmeta(y, S)
summary(mod1)

# CONSISTENCY MODEL, STRUCTURED BETWEEN-STUDY (CO)VARIANCE (PROPORTIONAL)
mod2 <- mixmeta(y, S, bscov="prop", control=list(Psifix=diag(3)+1))
summary(mod2)

# TRANSFORM IN LONG FORMAT, WITH S AS LIST (EXCLUDING MISSING)
long <- na.omit(reshape(smoking[,c(1,2,11:13)], varying=list(3:5), idvar="study", 
  v.names="y", timevar="outcome", times=colnames(y), direction="long"))
Slist <- lapply(lapply(seq(nrow(S)), function(i) xpndMat(S[i,])), function(x)
  x[!is.na(diag(x)), !is.na(diag(x)), drop=F])

# THE MODELS ABOVE CAN BE REPLICATED IN THE LONG FORMAT
mod2b <- mixmeta(y ~ 0 + factor(outcome), random= ~ 0 + factor(outcome)|study,
  data=long, bscov="prop", control=list(addS=Slist, Psifix=diag(3)+1))
summary(mod2b)

# DEFINE AND ADD INDICATORS FOR OUTCOME AND DESIGN
dummy <- cbind(model.matrix(~0+outcome, long), model.matrix(~0+design, long))
colnames(dummy) <- c(levels(factor(long$outcome)), levels(long$design))
long <- cbind(long, data.frame(dummy))

# INCONSISTENCY MODEL (SPECIAL PARAMETERIZATION OF OUTCOME-BY-DESIGN INTERACTION)
formula <- y ~ 0 + yB + yC + yC:acd + yC:bc + yC:bcd + yD + yD:acd + yD:bcd + 
  yD:bd + yD:cd
mod3 <- update(mod2b, formula=formula)
summary(mod3)

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}
fwald(mod3, c(":"))

# SEE help(smoking) FOR FURTHER INFO
