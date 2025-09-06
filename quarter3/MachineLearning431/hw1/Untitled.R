library(ggplot2)
library(olsrr)

### Question 1 ###
data(diamonds)
cutf=as.character(diamonds$cut)
cutf=as.factor(cutf)
qplot(carat, price, data = diamonds, col = cutf, size = I(0.5)) + theme_bw()
qplot(carat, log(price), data = diamonds, col = cutf, size = I(0.5)) + theme_bw()
qplot(carat, sqrt(price), data = diamonds, col = cutf, size = I(0.5)) + theme_bw()

model = lm(log(price) ~ carat, data = diamonds)
#ols_plot_resid_fit(model)
qplot(model$fitted.values, model$residuals, col = I("blue")) + geom_abline(slope = 0, intercept = 0, col = I("red"))
qplot(diamonds$carat, model$residuals, col = I("blue")) + geom_abline(slope = 0, intercept = 0, col = I("red"))

# Most diamonds' weights concentrate on the range from 0 to 2.5 carats. In this range,
# residuals of this regression are basically around 0 and symmetrical, indicating that
# this is a approximate linear relationship. But when carat is greater than 2.5, there 
# emerges many outliers, with residuals way less than 0, which means for the whole range,
# the linear relationship is not perfect.


### Question 2 ###
library(data.table)
library(foreign)
#setwd("/Users/leonard/Desktop/MFE-ML/hw1/")
StockRetAcct_DT <- as.data.table(read.dta('/Users/huanyu/Desktop/MachineLearning431/hw1/StockRetAcct_insample.dta'))
#Sort portfolios based on lnIssue each year
for(i in 1980:2014){
    StockRetAcct_DT[year==i,Ins_decile_yr:=cut(StockRetAcct_DT[year==i,]$lnIssue,breaks=quantile(StockRetAcct_DT[year==i,]$lnIssue, probs=c(0:10)/10, na.rm = TRUE),include.lowest = TRUE, labels = FALSE)]
}
# Equal-weighted, average return each year for each decile portfolio.
EW_Ins_MutualFund_yr = StockRetAcct_DT[,list(MeanRetYear = mean(lnAnnRet)), by = list(Ins_decile_yr, year)]
# EV, average return across years
EW_Ins_MutualFund = EW_Ins_MutualFund_yr[,list(MeanRet = mean(MeanRetYear)), by = Ins_decile_yr]
# Value-weighted
VW_Ins_MutualFund_yr = StockRetAcct_DT[,list(MeanRetYear = weighted.mean(lnAnnRet, MEwt)), by = list(Ins_decile_yr, year)]
VW_Ins_MutualFund = VW_Ins_MutualFund_yr[,list(MeanRet = mean(MeanRetYear)), by = Ins_decile_yr]
setkey(EW_Ins_MutualFund, Ins_decile_yr)
setkey(VW_Ins_MutualFund, Ins_decile_yr)
# Compare VW and EW portfolios
Comparision = merge(EW_Ins_MutualFund,VW_Ins_MutualFund,by="Ins_decile_yr")

qplot(Ins_decile_yr, MeanRet, data = VW_Ins_MutualFund, na.rm = TRUE, col = I("blue"), main = "VW Insurance bins vs. Returns") + geom_smooth(method = "lm", col=I("red"))
qplot(Ins_decile_yr, MeanRet, data = EW_Ins_MutualFund, na.rm = TRUE, col = I("blue"), main = "EW Insurance bins vs. Returns") + geom_smooth(method = "lm", col=I("red"))

# Create a new columne of transIns
StockRetAcct_DT[,transIns:= 0][Ins_decile_yr==1,transIns:= -1][Ins_decile_yr==10,transIns:= 1][is.na(Ins_decile_yr), transIns:= NA]

LM <- function(x){
    model <- lm(lnAnnRet ~ transIns, data = x)
    coef(model)[2]
}
Lamda_t = StockRetAcct_DT[!is.na(Ins_decile_yr),LM(.SD), by = year]
setnames(Lamda_t,c("year","lamda"))
# lamda of Fama-MacBeth 
lamda = mean(Lamda_t$lamda,na.rm = T)
lamda
# Long decile 1, short decile 10

### Question 3 ###

for(i in 1980:2014){
    StockRetAcct_DT[year==i,size_quintile_yr:=cut(StockRetAcct_DT[year==i,]$MEwt,breaks=quantile(StockRetAcct_DT[year==i,]$MEwt, probs=c(0:5)/5, na.rm = TRUE),include.lowest = TRUE, labels = FALSE)]
    StockRetAcct_DT[year==i,bm_quintile_yr:=cut(StockRetAcct_DT[year==i,]$lnBM,breaks=quantile(StockRetAcct_DT[year==i,]$lnBM, probs=c(0:5)/5, na.rm = TRUE),include.lowest = TRUE, labels = FALSE)]
}
# Equal-weighted, average return each year for each Size bin and each BM bin.
EW_BM_MutualFund_yr = StockRetAcct_DT[,list(MeanRetYear = mean(lnAnnRet)), by = list(bm_quintile_yr, size_quintile_yr, year)]
# EV, average return across years
EW_BM_MutualFund = EW_BM_MutualFund_yr[,list(MeanRet = mean(MeanRetYear)), by = list(bm_quintile_yr, size_quintile_yr)]

qplot(bm_quintile_yr, MeanRet, facets = size_quintile_yr~., data=EW_BM_MutualFund) + geom_smooth(method = "lm", col=I("red"))