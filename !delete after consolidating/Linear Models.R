# https://ucdavis-bioinformatics-training.github.io/2019-March-Bioinformatics-Prerequisites/thursday/linear_models.html

#******************************************************************************#
# Linear Models
#******************************************************************************#

# A linear model is a model for a continuous outcome Y of the form
# Y = β0 + β1X1 + β2X2 + ⋯ + βpXp + ϵ

# The covariates X can be:
# (i) a continuous variable (age, weight, temperature, etc.)
# (ii) Dummy variables coding a categorical covariate

# The β’s are unknown parameters to be estimated.
# The error term ϵ is assumed to be normally distributed with a variance that is
# constant across the range of the data.

# Models with all categorical covariates are referred to as ANOVA models and
# models with continuous covariates are referred to as linear regression models.
# These are all linear models, and R doesn’t distinguish between them.

#******************************************************************************#
# Linear Models in R
#******************************************************************************#

library(ggplot2)

dat <- read.csv("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2018-September-Bioinformatics-Prerequisites/master/friday/lm_example_data.csv")
head(dat)

# Fit a linear model using expression as the outcome and treatment as a 
# categorical covariate. In R model syntax, the outcome is on the left side, 
# with covariates (separated by +) following the ~
oneway.model <- stats::lm(formula = expression ~ treatment, 
                          data = dat)

#create scatterplot with fitted regression line
ggplot2::ggplot(dat, aes(x = treatment, y = expression)) +
  geom_point() +
  stat_smooth(method = "lm") 

# View the summary of the linear fit
# “Coefficients” refer to the β’s
# “Estimate” is the estimate of each coefficient
# “Std. Error” is the standard error of the estimate
# “t value” is the coefficient divided by its standard error
# “Pr(>|t|)” is the p-value for the coefficient
# The residual standard error is the estimate of the variance of ϵ
# Degrees of freedom is the sample size minus number of coefficients estimated
# R-squared is (roughly) the proportion of variance in the outcome explained by the model
# The F-statistic compares the fit of the model as a whole to the null model (with no covariates)
base::summary(oneway.model)

# View the estimated values of the co-efficients
stats::coef(oneway.model)

# What do the model coefficients mean?
# By default, R uses reference group coding or “treatment contrasts”. 
# For categorical covariates, the 1st level alphabetically (or 1st factor level)
# is treated as the reference group. The reference group doesn’t get its own 
# coefficient, it is represented by the intercept. 

# Coefficients for other groups are the difference from the reference:
# For our simple design:
# (Intercept) is the mean of expression for treatment A
# treatment B is the mean of expression for treatment B minus the mean for treatment A
# treatment C is the mean of expression for treatment C minus the mean for treatment A

mean(dat$expression[dat$treatment == "A"])
mean(dat$expression[dat$treatment == "B"]) - mean(dat$expression[dat$treatment == "A"])
mean(dat$expression[dat$treatment == "C"]) - mean(dat$expression[dat$treatment == "A"])
mean(dat$expression[dat$treatment == "D"]) - mean(dat$expression[dat$treatment == "A"])
mean(dat$expression[dat$treatment == "E"]) - mean(dat$expression[dat$treatment == "A"])

# If you don’t want reference group coding, then fit a model without an intercept:
no.intercept.model <- stats::lm(formula = expression ~ 0 + treatment, 
                                data = dat) # '0' means 'no intercept' here
coef(no.intercept.model)

# Without the intercept, the coefficients here estimate the mean in each level of treatment:
mean(dat$expression[dat$treatment == "A"])
mean(dat$expression[dat$treatment == "B"])
mean(dat$expression[dat$treatment == "C"])
mean(dat$expression[dat$treatment == "D"])
mean(dat$expression[dat$treatment == "E"])

#******************************************************************************#
# The Design Matrix
#******************************************************************************#

# For the RNASeq analysis programs limma and edgeR, the model is specified 
# through the design matrix.
 
# The design matrix X has one row for each observation and one column for each 
# model coefficient. The design matrix can be specified through the model.matrix
# function using the same syntax as for lm, just without a response:
   
# Design matrix for reference group coded model:
X <- model.matrix(object = ~treatment, data = dat)

#******************************************************************************#
# Adding More Covariates
#******************************************************************************#

# Batch Adjustment
# Suppose we want to adjust for batch differences in our model. We do this by 
# adding the covariate “batch” to the model formula:
  
batch.model <- lm(expression ~ treatment + batch, data = dat)
summary(batch.model)
batch.X <- model.matrix(object = ~treatment + batch, data = dat)

#******************************************************************************#
# Two-Way ANOVA Models
#******************************************************************************#

# Suppose our experiment involves two factors, treatment and time. lm can be 
# used to fit a two-way ANOVA model:

twoway.model <- lm(expression ~ treatment*time, data = dat)
summary(twoway.model)
twoway.X <- model.matrix(object = ~treatment*time, data = dat)

# The notation treatment*time refers to treatment, time, and the interaction 
# effect of treatment by time.(This is different from other statistic software).

# Interpretation of coefficients:
# Each coefficient for treatment represents the difference between the indicated
# group and the reference group at the reference level for the other covariates
# For example, 
# “treatmentB” = difference in expression between treatment B and treatment A at time 1
# “timetime2”  = difference in expression between time2 and time1 for treatment A
# The interaction effects (coefficients with “:”) estimate 
# (i) the difference between treatment groups in the effect of time
# (ii) the difference between times in the effect of treatment

# Another Parameterization
# In a multifactor model, estimating contrasts can be fiddly, especially with 
# lots of factors or levels. Here is an equivalent way to estimate the same 
# two-way ANOVA model that gives easier contrasts:
  
# First, define a new variable that combines the information from treatment and 
# time variables

dat$tx.time <- base::interaction(dat$treatment, dat$time)
dat

# Next, fit a one-way ANOVA model with the new covariate. Don’t include an 
# intercept in the model.

other.2way.model <- lm(expression ~ 0 + tx.time, data = dat)
summary(other.2way.model)

# We get the same estimates from both two way models
# Effect of treatment B vs. A at time 1:
c1 <- coef(twoway.model)
c1["treatmentB"] 

c2 <- coef(other.2way.model)
c2["tx.timeB.time1"] - c2["tx.timeA.time1"]

# Effect of treatment B vs. A at time 2:
c1 <- coef(twoway.model)
c1["treatmentB"] + c1["treatmentB:timetime2"]

c2 <- coef(other.2way.model)
c2["tx.timeB.time2"] - c2["tx.timeA.time2"]

# Interaction effect (remembering that an interaction effect here is a difference of differences):
  
c1 <- coef(twoway.model)
c1["treatmentB:timetime2"]

c2 <- coef(other.2way.model)
(c2["tx.timeB.time2"] - c2["tx.timeA.time2"]) - (c2["tx.timeB.time1"] - c2["tx.timeA.time1"])

#******************************************************************************#
# Continuous Covariates
#******************************************************************************#

# Linear models with continuous covariates (“regression models”) are fitted in 
# much the same way:
  
continuous.model <- lm(expression ~ temperature, data = dat)
summary(continuous.model)
coefs <- coef(continuous.model)

# For the above model, the intercept is the expression at temperature 0 and the
# “temperature” coefficient is the slope, or how much expression increases for 
# each unit increase in temperature:

# Create scatterplot with fitted regression line
ggplot(dat, aes(x = temperature, y = expression)) +
  geom_point() +
  stat_smooth(method = "lm") +
  annotate("text", x=12, y=10, label= paste0("expression = ", round(coefs[1], 2),  "+", round(coefs[2], 2), "*temperature"))

# The slope from a linear regression model is related to but not identical to
# the Pearson correlation coefficient
cor.test(dat$expression, dat$temperature)

# Notice that the p-values for the correlation and the regression slope are 
# identical.

# Scaling and centering both variables yields a regression slope equal to the 
# correlation coefficient:
scaled.mod <- lm(scale(expression) ~ scale(temperature), data = dat)
coef(scaled.mod)[2]
