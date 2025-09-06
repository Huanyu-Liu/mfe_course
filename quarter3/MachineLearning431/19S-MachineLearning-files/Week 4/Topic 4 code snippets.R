#
# code snippets for Topic 4: Model Selection and Shrinkage
#

setwd("D:/llochsto/Dropbox/Data Analytics/Data")

require(glmnet)
require(foreign)
require(data.table)

# read in and prepare the Lending Club Data
LC_RawData <- as.data.table(read.dta("D:/llochsto/Dropbox/Data Analytics/Data/LendingClub_LoanStats3a_v12.dta"))
LC_Baseline <- LC_RawData[(loan_status == "Fully Paid") | (loan_status == "Charged Off")]
head(LC_Baseline)
LC_Baseline <- LC_Baseline[,Default:=(loan_status=="Charged Off")]

# create function that makes dummies out of factor in data.table
# we need to do this as glmnet does not automatically implement this in the same way that glm does
factorToDummy <- function(dtable, var.name){
  stopifnot(is.data.table(dtable))
  stopifnot(var.name %in% names(dtable))
  stopifnot(is.factor(dtable[, get(var.name)]))
  
  dtable[, paste0(var.name,"_",levels(get(var.name)))] -> new.names
  dtable[, (new.names) := transpose(lapply(get(var.name), FUN = function(x){x == levels(get(var.name))})) ]
  
  cat(paste("\nDummies created: ", paste0(new.names, collapse = ", ")))
}

# make sure variable of interest is a factor
LC_Baseline <- LC_Baseline[,grade_factor:=as.factor(grade)]

# create dummies, needed for glmnet
factorToDummy(LC_Baseline,"grade_factor")
head(LC_Baseline)

# create regressor matrix (glmnet needs matrix as x-input and vector as y-input)
regressors <- as.matrix(LC_Baseline[,list(grade_factor_A,grade_factor_B,grade_factor_C,grade_factor_D,grade_factor_E,grade_factor_F,grade_factor_G)])

# run the logistic ridge regression. Ridge is done when alpha = 0.
outloans_ridge = cv.glmnet(regressors,as.vector(LC_Baseline$Default), family=c("binomial"),alpha = 0, standardize = TRUE)

# plot the coefficients, as a function of the log of the constraint value given as "lambda".
# Note, this is Lagrangian formulation, high value of constraint means it binds more
plot.glmnet(outloans_ridge$glmnet.fit, "lambda", label = TRUE)
plot.cv.glmnet(outloans_ridge)


# non-grade model
LC_Baseline <- LC_Baseline[,term_factor:=as.factor(term)]
factorToDummy(LC_Baseline,"term_factor")
regressors <- as.matrix(LC_Baseline[,list(loan_amnt/10000,`term_factor_ 36 months`, `term_factor_ 60 months`, annual_inc/1000, 10*int_rate)])
outloans_ridge = cv.glmnet(regressors,as.vector(LC_Baseline$Default), family=c("binomial"),alpha = 0, standardize = TRUE)
plot.glmnet(outloans_ridge$glmnet.fit, "lambda", label = TRUE)
plot.cv.glmnet(outloans_ridge)


# including both grade factors and existing
regressors <- cbind(regressors, as.matrix(LC_Baseline[,list(grade_factor_A,grade_factor_B,grade_factor_C,grade_factor_D,grade_factor_E,grade_factor_F,grade_factor_G)]))
outloans_ridge = cv.glmnet(regressors,as.vector(LC_Baseline$Default), family=c("binomial"),alpha = 0, standardize = TRUE)
plot.glmnet(outloans_ridge$glmnet.fit, "lambda", label = TRUE)
plot.cv.glmnet(outloans_ridge)


# predict for new person with set of x attributes
new_regs <- matrix(c(50000/10000,1,0,200000/10000,0.05*10,0,1,0,0,0,0,0),nrow = 1)
predict(outloans_ridge, new_regs, s="lambda.1se", type = "response")



# now, let's turn to the lasso
outloans_lasso = cv.glmnet(regressors,as.vector(LC_Baseline$Default), family=c("binomial"),alpha = 1, standardize = TRUE)
plot.glmnet(outloans_lasso$glmnet.fit, "lambda", label = TRUE)
plot.cv.glmnet(outloans_lasso)
predict(outloans_lasso, new_regs, s="lambda.1se", type = "response")



