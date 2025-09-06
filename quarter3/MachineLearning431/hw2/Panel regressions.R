#
# code snippets for Topic 2: Panel regressions
#

# set your working directory, where you have downloaded the Stata data file
# change the below to match your folder
setwd("D:/llochsto/Dropbox/Data Analytics/DAML_2018/Data")


# we need the foreign package to import data in different format
require(foreign)
require(data.table)
require(ggplot2)

# Download data and set as data.table
StockRetAcct_DT <- as.data.table(read.dta("StockRetAcct_insample.dta"))

# set keys for StockRetAcct_DT In particular, it will be useful to sort on FirmID and year
setkey(StockRetAcct_DT, FirmID, year)

# create excess returns (what we really care about)
StockRetAcct_DT[,ExRet:=exp(lnAnnRet) - exp(lnRf)]


# Fama-MacBeth Regressions
# loop through the years in the data base
port_ret = NULL
for (i in 1980:2014)
{
  temp <- StockRetAcct_DT[year==i,]
  fit_yr <- lm(temp$ExRet ~ temp$lnBM, data=temp)
  temp <- coefficients(fit_yr)
  port_ret = rbind(port_ret,temp[2])
}
fm_output = list(MeanReturn = mean(port_ret), StdReturn = sqrt(var(port_ret)), SR_Return = mean(port_ret)/sqrt(var(port_ret)), tstat_MeanRet = sqrt(1+2014-1980)*mean(port_ret)/sqrt(var(port_ret)))
fm_output


# Let's try instead with another anomaly -- profitability, lnProf
port_ret = NULL
for (i in 1980:2014)
{
  temp <- StockRetAcct_DT[year==i,]
  fit_yr <- lm(temp$ExRet ~ temp$lnProf, data=temp)
  temp <- coefficients(fit_yr)
  port_ret = rbind(port_ret,temp[2])
}
fm_output = list(MeanReturn = mean(port_ret), StdReturn = sqrt(var(port_ret)), SR_Return = mean(port_ret)/sqrt(var(port_ret)), tstat_MeanRet = sqrt(1+2014-1980)*mean(port_ret)/sqrt(var(port_ret)))
fm_output
# profitability-based sort is significant, however... 

# What do you think the relation  between book-to-market and profitability is?
prof_bm = NULL
for (i in 1980:2014)
{
  temp <- StockRetAcct_DT[year==i,]
  fit_yr <- lm(temp$lnProf ~ temp$lnBM, data=temp)
  temp <- coefficients(fit_yr)
  prof_bm = rbind(prof_bm,temp[2])
}
fm_output = list(MeanCoeff = mean(prof_bm), tstat_Coeff = sqrt(1+2014-1980)*mean(prof_bm)/sqrt(var(prof_bm)))
fm_output
# highly significantly negative relation. Makes sense. High profitability should mean low book-to-market


# let's try both lnBM and lnProf together
port_ret = NULL
for (i in 1980:2014)
{
  temp <- StockRetAcct_DT[year==i,]
  fit_yr <- lm(temp$ExRet ~ temp$lnBM + temp$lnProf, data=temp)
  temp <- coefficients(fit_yr)
  port_ret = rbind(port_ret,temp[2:length(temp)])
}
fm_output = list(MeanReturn = colMeans(port_ret), StdReturn = sqrt(diag(var(port_ret))), SR_Return = colMeans(port_ret)/sqrt(diag(var(port_ret))), tstat_MeanRet = sqrt(1+2014-1980)*colMeans(port_ret)/sqrt(diag(var(port_ret))))
fm_output
# notice how in the multiple regression lnBM is slightly more significant, but with a lower mean. Why did this change?

# let's add industry dummies
port_ret = NULL
for (i in 1980:2014)
{
  temp <- StockRetAcct_DT[year==i,]
  #temp <- temp(,ind_fact)
  fit_yr <- lm(temp$ExRet ~ temp$lnProf + temp$lnBM + as.factor(temp$ff_ind), data=temp)
  temp <- coefficients(fit_yr)
  port_ret = rbind(port_ret,temp[2:length(temp)])
}
fm_output = list(MeanReturn = colMeans(port_ret), StdReturn = sqrt(diag(var(port_ret))), SR_Return = colMeans(port_ret)/sqrt(diag(var(port_ret))), tstat_MeanRet = sqrt(1+2014-1980)*colMeans(port_ret)/sqrt(diag(var(port_ret))))
fm_output
# now, both profitability and value are signficiant. What happened?


# getting portfolio weights for "new" value factor implied by the preceding FMB regressions
# choose the current date (end of sample) 
LastDate <- na.omit(StockRetAcct_DT[year==2014,])
Nt = length(LastDate$ExRet)
Xmat <- cbind(rep(1,Nt),LastDate$lnProf,LastDate$lnBM)
for (ii in 1:11) # note drop last industry dummy as we have intercept in regression
{
  Xmat <- cbind(Xmat,(LastDate$ff_ind==ii)) 
}
portweights_lnBM = solve(t(Xmat)%*%Xmat)%*%t(Xmat)
# lnBM is third row (first is intercept, second is profitability, fourth and on are industry)
portweights_lnBM = c(0,0,1,0,0,0,0,0,0,0,0,0,0,0) %*% portweights_lnBM
# scale portfolio weights to get 15% standard deviation of returns
portweights_lnBM = portweights_lnBM * 0.15 / sqrt(var(port_ret[,2]))

# for plotting, get the scaled excess portfolio returns
lnBM_ret = port_ret[,2] * 0.15 / sqrt(var(port_ret[,2]))
# create cumulative log return series
cum_ret_lnBM = 0
for (ii in 1:35)
{
  cum_ret_lnBM = rbind(cum_ret_lnBM,cum_ret_lnBM[ii]+log(1+lnBM_ret[ii]))
}

# get "old" simple value strategy returns
port_ret = NULL
for (i in 1980:2014)
{
  temp <- StockRetAcct_DT[year==i,]
  fit_yr <- lm(temp$ExRet ~ temp$lnBM, data=temp)
  temp <- coefficients(fit_yr)
  port_ret = rbind(port_ret,temp[2])
}
lnBM_old_ret = port_ret[,1] * 0.15 / sqrt(var(port_ret[,1]))
cum_ret_oldlnBM = 0
for (ii in 1:35)
{
  cum_ret_oldlnBM = rbind(cum_ret_oldlnBM,cum_ret_oldlnBM[ii]+log(1+lnBM_old_ret[ii]))
}
# plot exponential of cumulkative log return to get to regular cumulative returns
# shows pretty convincly how the new cleaned-up value strategy performs better than the old
qplot(c(1980:2015), exp(cum_ret_oldlnBM), geom="line", xlab="year",ylab="Cumulative Return",color = I("blue"), size=I(1.5), main = "Old Value (blue) vs. New Value (red)") +
  geom_line(aes(y=exp(cum_ret_lnBM)),color = I("red"),size=I(1.5)) +
  theme_bw() 




# EARNINGS FORECASTING IN PANEL REGRESSIONS
# Prep
library(foreign)
library(data.table)
library(lfe)
library(stargazer)
rm(list = ls())

# Read data
StockRetAcct_DT = data.table(read.dta("StockRetAcct_insample.dta"))

# create excess returns (what we really care about)
StockRetAcct_DT[,ExRet:=exp(lnAnnRet) - exp(lnRf)]

# Get mean and standard deviation of ROE across firms and time
mean_ROE = mean(na.omit(StockRetAcct_DT$lnROE))
mean_ROE
std_ROE = sqrt(var(na.omit(StockRetAcct_DT$lnROE)))
std_ROE

# What predicts next year's ROE?
setorder(StockRetAcct_DT, FirmID, year) # Set order of data so that shift does what we want it to
StockRetAcct_DT[, lead_lnROE := shift(lnROE, type = 'lead'), by = FirmID]  # Define next year's lnROE
roe_panel = felm(lead_lnROE ~ lnROE, StockRetAcct_DT) # Regression with no FE or SE
stargazer(roe_panel, type = 'text', report = 'vc*t') # Output regressions as text, report t-stats

# cluster standard errors at the firm level
roe_panel2 = felm(lead_lnROE ~ lnROE | 0 | 0 | FirmID, StockRetAcct_DT) # Regression with no FE, standard errors clustered at the firm level
stargazer(roe_panel2, type = 'text', report = 'vc*t') # Output regressions as text, report t-stats

# clustered standard errors at the firm and year level
roe_panel3 = felm(lead_lnROE ~ lnROE | 0 | 0 | year + FirmID, StockRetAcct_DT) # Regression with no FE, standard errors clustered at the firm and year level
stargazer(roe_panel3, type = 'text', report = 'vc*t') # Output regressions as text, report t-stats

# show output in one table, neat feature of stargazer
stargazer(roe_panel, roe_panel2, roe_panel3, type = 'text', report = 'vc*t') # Output regressions as text, report t-stats

# industry component in accounting variables, include industry fixed effects
roe_panel4 = felm(lead_lnROE ~ lnROE| ff_ind | 0 | year + FirmID, StockRetAcct_DT) # Regression with FE, standard errors clustered at the firm and year level
stargazer(roe_panel4, type = 'text', report = 'vc*t') # Output regressions as text, report t-stats

# industry component in accounting variables, include industry and year fixed effects
roe_panel5 = felm(lead_lnROE ~ lnROE| year + ff_ind | 0 | year + FirmID, StockRetAcct_DT) 
stargazer(roe_panel4, roe_panel5, type = 'text', report = 'vc*t') # Output regressions as text, report t-stats


# clustered standard errors at the firm and year level, industry effects: big model
roe_panel6 = felm(lead_lnROE ~ lnROE + lnBM +lnProf + lnLever + lnIssue + lnInv | ff_ind | 0 | year + FirmID, StockRetAcct_DT)
stargazer(roe_panel6, type = 'text', report = 'vc*t') # Output regressions as text, report t-stats


# What predicts ROE in five  years? 
StockRetAcct_DT[, lead2_lnROE := shift(lead_lnROE, type = 'lead'), by = FirmID]  # Define next year's lnROE
StockRetAcct_DT[, lead3_lnROE := shift(lead2_lnROE, type = 'lead'), by = FirmID]  # Define next year's lnROE
StockRetAcct_DT[, lead4_lnROE := shift(lead3_lnROE, type = 'lead'), by = FirmID]  # Define next year's lnROE
StockRetAcct_DT[, lead5_lnROE := shift(lead4_lnROE, type = 'lead'), by = FirmID]  # Define next year's lnROE
roe_panel7 = felm(lead5_lnROE ~ lnROE + lnBM +lnProf + lnLever + lnIssue + lnInv | ff_ind | 0 | year + FirmID, StockRetAcct_DT) # Regression with no FE, standard errors clustered at the firm level
stargazer(roe_panel6, roe_panel7, type = 'text', report = 'vc*t') # Output regressions as text, report t-stats



