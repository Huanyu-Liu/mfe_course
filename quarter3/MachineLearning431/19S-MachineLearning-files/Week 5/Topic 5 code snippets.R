#
# code snippets for Topic 5: Decision Trees
# code shows example of RandomForest on the usual stock return data set
#

setwd("C:/Users/llochsto/Dropbox/Data Analytics/Data")
#setwd("D:/llochsto/Dropbox/Data Analytics/Data")


require(glmnet)
require(foreign)
require(data.table)
require(randomForest)

# Download data and set as data.table
StockRetAcct_DT <- as.data.table(read.dta("StockRetAcct_insample.dta"))

# take a look at content, first 6 rows. Notice there are missings, data is annual.
head(StockRetAcct_DT)

# set keys for StockRetAcct_DT In particular, it will be useful to sort on FirmID and year
setkey(StockRetAcct_DT, FirmID, year)
setorder(StockRetAcct_DT, year, FirmID)


# create excess returns (what we really care about)
StockRetAcct_DT[,ExRet:=exp(lnAnnRet) - exp(lnRf)]

# remove any rows that contain missing data
StockRetAcct_DT <- StockRetAcct_DT[complete.cases(StockRetAcct_DT),]

# train model on data up until 2010
y_data_train <- as.vector(StockRetAcct_DT[year<2010,ExRet])
x_data_train <- as.matrix(StockRetAcct_DT[year<2010,list(lnIssue,lnMom,lnME,lnProf,lnEP,lnInv,lnLever,lnROE,rv,lnBM,ff_ind,year)])

# test data is data from 2010 and on 
y_data_test <- as.vector(StockRetAcct_DT[year>=2010,ExRet])
x_data_test <- as.matrix(StockRetAcct_DT[year>=2010,list(lnIssue,lnMom,lnME,lnProf,lnEP,lnInv,lnLever,lnROE,rv,lnBM,ff_ind,year)])

# set maximum terminal nodes equal to 30 for a relatively parsimonious trees
rfLC <- randomForest(x_data_train, y_data_train, ntree = 500, maxnodes = 30)

# check prediction fit in- and out-of-sample
RF_pred_in_sample <- predict(rfLC, x_data_train)
RF_test_in_sample <- lm(y_data_train~RF_pred_in_sample)
summary(RF_test_in_sample)
RF_pred_out_of_sample <- predict(rfLC, x_data_test)
RF_test_out_of_sample <- lm(y_data_test~RF_pred_out_of_sample)
summary(RF_test_out_of_sample)

# compare to linear regression
x_lm <- x_data_train
LR_test_in_sample <- lm(y_data_train~x_lm)
summary(LR_test_in_sample)

x_lm <- x_data_test
LR_pred_out_of_sample <- predict(LR_test_in_sample, as.data.frame(x_lm))
LR_test_out_of_sample <- lm(y_data_test~LR_pred_out_of_sample)
summary(LR_test_out_of_sample)

# Fama-MacBeth regressions using out-of-sample predicted values
# Random Forest first
port_ret = NULL
row_counter_end = 0
for (i in 2010:2014)
{
  year_length <- StockRetAcct_DT[year==i,ExRet]
  year_length <- length(year_length)
  row_counter_start = row_counter_end + 1
  row_counter_end = row_counter_end + year_length
  x_temp <- RF_pred_out_of_sample[row_counter_start:row_counter_end]
  y_temp <- y_data_test[row_counter_start:row_counter_end]
  fit_yr <- lm(y_temp ~ x_temp)
  summary(fit_yr)
  temp <- coefficients(fit_yr)
  port_ret = rbind(port_ret,temp[2])
}
# recall, scale depends on magnitude of x-variable, so only fair to compare Sharpe ratios and t-stats
fm_RF_output = list(SR_Return = mean(port_ret)/sqrt(var(port_ret)), tstat_MeanRet = sqrt(1+2014-2010)*mean(port_ret)/sqrt(var(port_ret)))
fm_RF_output


# Linear regression next
port_ret = NULL
row_counter_end = 0
for (i in 2010:2014)
{
  year_length <- StockRetAcct_DT[year==i,ExRet]
  year_length <- length(year_length)
  row_counter_start = row_counter_end + 1
  row_counter_end = row_counter_end + year_length
  x_temp <- LR_pred_out_of_sample[row_counter_start:row_counter_end]
  y_temp <- y_data_test[row_counter_start:row_counter_end]
  fit_yr <- lm(y_temp ~ x_temp)
  temp <- coefficients(fit_yr)
  port_ret = rbind(port_ret,temp[2])
}
fm_LR_output = list(SR_Return = mean(port_ret)/sqrt(var(port_ret)), tstat_MeanRet = sqrt(1+2014-2010)*mean(port_ret)/sqrt(var(port_ret)))
fm_LR_output


# Next, let's see how the XGBoost method does in the same exercise
require(xgboost)

# first, let's set parameters for the xgboost, no regularization (gamma = 0) other than the cv procedure
params <- list(booster = "gbtree", objective = "reg:linear", eta = 0.3, gamma = 0, max_depth = 20)

# XGBoost likes to use xgb.DMatrix
xgb_train <- xgb.DMatrix(data = x_data_train, label = y_data_train)

# use xgb.cv to find the best nround (number of trees) for this model.
xgbcv <- xgb.cv(params = params, data = xgb_train, nfold = 10, nrounds = 100, showsd = T, print_every_n = 10) 

# pick out the lowest cv mean rmse
cv_nrounds = which.min(xgbcv$evaluation_log$test_rmse_mean)

# with optimal nrounds in hand, run the prediction for out of sample
xgb_optb <- xgboost(params = params, data = xgb_train, nround = cv_nrounds)
xgb_test <- xgb.DMatrix(data = x_data_test, label = y_data_test)
xgb_pred <- predict(xgb_optb, xgb_test)
xgb_test_out_of_sample <- lm(y_data_test~xgb_pred)
summary(xgb_test_out_of_sample)

# plot importance of features in order to understand model
boost_out <- xgb.importance(feature_names = colnames(x_data_test), model = xgb_optb)
xgb.plot.importance(importance_matrix = boost_out[1:11])


# FMB regression to get tradeable portfolio
# Linear regression next
port_ret = NULL
row_counter_end = 0
for (i in 2010:2014)
{
  year_length <- StockRetAcct_DT[year==i,ExRet]
  year_length <- length(year_length)
  row_counter_start = row_counter_end + 1
  row_counter_end = row_counter_end + year_length
  x_temp <- xgb_pred[row_counter_start:row_counter_end]
  y_temp <- y_data_test[row_counter_start:row_counter_end]
  fit_yr <- lm(y_temp ~ x_temp)
  temp <- coefficients(fit_yr)
  port_ret = rbind(port_ret,temp[2])
}
fm_LR_output = list(SR_Return = mean(port_ret)/sqrt(var(port_ret)), tstat_MeanRet = sqrt(1+2014-2010)*mean(port_ret)/sqrt(var(port_ret)))
fm_LR_output
