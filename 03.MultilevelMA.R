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
# MULTILEVEL META-ANALYSIS (SECTION 4.3 - EXAMPLE 1, TABLE 2)
################################################################################

# LOAD THE PACKAGE
library(mixmeta)

# STUDY AS SINGLE LEVEL: STANDARD META-ANALYSIS
mod1 <- mixmeta(effect, var, random= ~ 1|study, data=school, method="ml")
print(summary(mod1), digits=3, report="var")

# DISTRICT AS SINGLE LEVEL: META-ANAYSIS WITH REPEATED MEASURES
mod2 <- mixmeta(effect, var, random= ~ 1|district, data=school, method="ml")
print(summary(mod2), digits=3, report="var")

# NESTED LEVELS OF STUDY AND DISTRICT: TWO-LEVEL META-ANALYSIS
mod3 <- mixmeta(effect, var, random= ~ 1|district/study, data=school,
  method="ml")
print(summary(mod3), digits=3, report="var")

# COMPARISON
AIC(mod1, mod2, mod3)

# SEE help(school) FOR FURTHER EXAMPLES


################################################################################
# MULTILEVEL META-ANALYSIS (SECTION 4.3 - EXAMPLE 2, TABLE 3)
################################################################################

# STANDARD META-ANALYSIS: IGNORING CLUSTERING OF TRIALS
subtrial <- seq(nrow(thrombolytic))
mod1 <- mixmeta(absrisk, var, random= ~ 1|subtrial, data=thrombolytic)
print(summary(mod1), digits=5)

# STANDARD META-REGRESSION
mod2 <- mixmeta(absrisk~time2treat, var, random= ~ 1|subtrial,
  data=thrombolytic)
print(summary(mod2), digits=5)

# TWO-LEVEL META-ANALYSIS
mod3 <- mixmeta(absrisk, var, random= ~ 1|trial/subtrial, data=thrombolytic)
print(summary(mod3), digits=5)

# TWO-LEVEL META-REGRESSION 
mod4 <- mixmeta(absrisk~time2treat, var, random= ~ 1|trial/subtrial,
  data=thrombolytic)
print(summary(mod4), digits=5)

# SEE help(thrombolytic) FOR FURTHER INFO
