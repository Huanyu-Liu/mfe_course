rm(list = ls())
### Question 1 ###
library(ggplot2)
library(data.table)
library(readxl)
library(reshape2)
library(zoo)
library(randomForest)
library(lfe)
library(xgboost)

# Part A
dtable = as.data.table(read_xlsx("/Users/huanyu/Desktop/MachineLearning431/hw5/French_Portfolio_Returns.xlsx",sheet='ExcessPortRet',skip=1))
dtable = dtable[Date>=196307,]
dtable$Month = c(1:length(dtable$Date))
dtable = dtable[,-'Date']
dtable = as.data.frame(dtable) # transform into df to apply melt
dtable_cv = melt(data = dtable, id.vars = 'Month', variable.name = 'Portfolio',value.name = 'ExRet',na.rm = TRUE)
dtable_cv$Portfolio = as.numeric(dtable_cv$Portfolio)
dtable_cv$ExRet = as.numeric(dtable_cv$ExRet)
dtable_cv = as.data.table(dtable_cv)
dtable_cv[,lag_1:=shift(ExRet), by = Portfolio]
dtable_cv[,lag_2:=shift(ExRet,n = 2), by = Portfolio]
dtable_cv[,lag_3_12:=shift(rollapplyr(ExRet, 10, sum, fill='NA'),3), by = Portfolio]
dtable_cv = na.omit(dtable_cv)
dtable_cv[,lag_1_square:=lag_1**2]
dtable_cv[,lag_2_square:=lag_2**2]
dtable_cv[,lag_3_12_square:=lag_3_12**2]

# Part B
dat_train = dtable_cv[Month < 559]
y_data_train <- as.vector(dat_train$ExRet)
x_data_train <- as.matrix(dat_train[,c(4:9)])
rfLC <- randomForest(x_data_train, y_data_train, ntree = 500, maxnodes = 30)
RF_pred_in_sample <- predict(rfLC, x_data_train)
RF_test_in_sample = felm(y_data_train ~ RF_pred_in_sample)
summary(RF_test_in_sample)
RF_panel_in_sample = felm(ExRet ~ RF_pred_in_sample | 0 | 0 | Month, dat_train)
summary(RF_panel_in_sample)

# test data is data from 2010 and on 
dat_test = dtable_cv[Month > 558]
y_data_test <- as.vector(dat_test$ExRet)
x_data_test <- as.matrix(dat_test[,c(4:9)])
RF_pred_out_of_sample <- predict(rfLC, x_data_test)
RF_test_out_of_sample <- lm(y_data_test~RF_pred_out_of_sample)
summary(RF_test_out_of_sample)
RF_panel_out_of_sample = felm(ExRet ~ RF_pred_out_of_sample | 0 | 0 | Month, dat_test)
summary(RF_panel_out_of_sample)

# Fama-MacBeth regressions using out-of-sample predicted values
# Random Forest first
port_ret = NULL
row_counter_end = 0
for (i in 559:631){
    year_length <- dtable_cv[Month==i,ExRet]
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
fm_RF_output = list(SR_Return = mean(port_ret)/sqrt(var(port_ret)), tstat_MeanRet = sqrt(1+631-559)*mean(port_ret)/sqrt(var(port_ret)))
fm_RF_output
