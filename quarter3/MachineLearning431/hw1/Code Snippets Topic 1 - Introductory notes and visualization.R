#
# code snippets for Topic 1: Introductory nots and visualization
#
# key commands in base R graphics
#
par(mfrow=c(1,2))    # create a 1x2 array of plots, paint row by row from upper left
x=rnorm(1000)
y=rnorm(1000)+2
hist(x,main="Histogram of x",breaks=20,col="magenta",xlab="")
plot(x,y,xlab="x", ylab="y",pch=20,col="red",xlim=c(-4,4),ylim=c(-4,4))
abline(c(0,1),lwd=2,lty=2)
# points(z,w)  add points with coordinates z,w to the plot
# lines(z,w,lwd=2,lty=2) draw lines between points spec by z,w of width 2 and type 2 (dashed)
#
# try demo
#
demo(graphics)
#
# now let's try some features of ggplot2
library(ggplot2)
# make sure you also have installed the colorspace-package
#
# some simple plots

# first, install DataAnalytics package
install.packages("devtools")                            # install devtools package
library("devtools")                                     # "load" package
install_bitbucket("perossichi/DataAnalytics")           # install the course dataset package
# note: this will take a few seconds
# depending on internet connection speed
library("DataAnalytics")

#
data(mvehicles)
cars=mvehicles[mvehicles$bodytype != "Truck",]
qplot(log(emv),data=cars)
qplot(luxury,log(emv),data=cars)
#
# now let's put some different labels on the axis and a title and make background white and show
# some other possible options
#
qplot(luxury,log(emv),data=cars,ylab="log(Market Value)",col=I("blue"),
      main="log(Market Value) Vs. Luxury") + theme_bw()
# qplot is a function which produces plot meta-data and theme_bw() adds or modifies this data
#   behind the scenes there is a "rendering" function the paints the plot.
# xlab,ylab change labels on axes
# main sets the title
# col or colour associates col with a variable. 
#   If you want just a constant color for all points,
#   you need I("<color name>")

#
# let's try to visualize another variable, bodytype, in the same plot
#
qplot(luxury,log(emv),data=cars,col=bodytype) + theme_bw()
# or make the points bigger and "transparent"
qplot(luxury,log(emv),data=cars,col=bodytype,size=I(4),alpha=I(3/4)) + theme_bw()
#
# let's make a separate plot for each value in bodytype
qplot(luxury,log(emv),data=cars,facets=bodytype~.,col=I("blue")) + theme_bw()
#
# let's add a fitted line to each plot
#
qplot(luxury,log(emv),data=cars,facets=bodytype~.,col=I("blue"), 
      geom=c("point","smooth"), method="lm") + theme_bw()
#
# optional: I want the fitted lines to be thicker and in red
#           and I don't like the grey error band
#
qplot(luxury,log(emv),data=cars,facets=bodytype~.,col=I("blue")) +
  geom_smooth(method="lm",col=I("red"),size=I(1.2),se=FALSE) + theme_bw()

#
# relate a continuous variable like emv to a categorical variable like bodytype
#
qplot(bodytype,log(emv),data=cars,geom="boxplot",fill=I("green")) + theme_bw()



# data.table code snippets
# now let's learn about the data.table package
#
library(data.table)
sales.table = data.table(module=c("a","a","c","c","c","b"),upc=c(2,3,4,5,5,1),
                         sales=round(runif(6)*10000,0))
sales.table
setkey(sales.table,module,upc)  # makes module and upc keys and "sorts" in place
# show all data.tables in the workspace
#
tables()

#
# accessing rows of data table -- pretty much the same
sales.table[1:3,]
sales.table[upc > 3,]    
# note you don't have to say "sales.table$upc" in the logical expression

# 
# now we can use the keys to access rows as well
#
sales.table["a"]
sales.table[.(c("a","c"))]

# 
# to refer to variables we can use the standard $ notation
#
sales.table$module
#
# or use the name of the variable in the "j" position
sales.table[,list(module,upc)]

#
# let's create a new variable which is the log of sales
#
sales.table[,lnsales:=log(sales)]
sales.table

#
# let's find total sales for each module.
#
sales.module=sales.table[,sum(sales),by=module]
setnames(sales.module,"V1","agg_sales")
sales.module

#
# now let's use a look-up table to add module descriptions
#
#
# a "lookup" table of module descriptions
#
module.description = data.table(module=c("a","b"),
                                description=c("meats","milk"))
setkey(module.description,module)
sales.module[module.description]

# or 
merge(sales.module,module.description,by="module")

# or

merge(sales.module,module.description)



# now, at sub-bullets e) and f) in Topic 1
# reading Stata database with stock returns and one-year lagged accounting variables. StockRetAcct.dta
# For variable descriptions, see StockRetAcct_README.txt

# set your working directory, where you have downloaded the Stata data file
# change the below to match your folder
setwd("D:/llochsto/Dropbox/Data Analytics/DAML_2018/Data")

# we need the foreign package to import data in different format
require(foreign)
require(data.table)
require(ggplot2)

# Download data and set as data.table
StockRetAcct_DT <- as.data.table(read.dta("StockRetAcct_insample.dta"))

# take a look at content, first 6 rows. Notice there are missings, data is annual.
head(StockRetAcct_DT)

# set keys for StockRetAcct_DT In particular, it will be useful to sort on FirmID and year
setkey(StockRetAcct_DT, FirmID, year)

# create excess returns (what we really care about)
StockRetAcct_DT[,ExRet:=exp(lnAnnRet) - exp(lnRf)]

# scatter plot excess returns vs. lagged bm 
# (notice that all firm characteristics are already lagged so we don't need to deal with this)
qplot(lnBM, ExRet,data = StockRetAcct_DT, col=I("blue"), main = "Log Book-to-Market vs. Excess Returns") + theme_bw()

# The data is messy. A slight positive relationship is arguably visible. 
# Let's fit a curve to this. Exclude missing observations for this to avoid error message:
qplot(lnBM, ExRet,data = StockRetAcct_DT[!is.na(StockRetAcct_DT$lnBM)], col=I("blue"), main = "Log Book-to-Market vs. Excess Returns") + geom_smooth(col=I("red"),se=FALSE,na.rm=TRUE) + theme_bw()

# the fit is done with additive generalized linear models. 
# Let's zoom in without removing data points (i.e, don't use ylim and xlim) and make the fitting function use more local data
# thus, allowing for more nonlinear fit by using the "span" option
qplot(lnBM, ExRet,data = StockRetAcct_DT[!is.na(StockRetAcct_DT$lnBM)], col=I("blue"), main = "Log Book-to-Market vs. Excess Returns") + geom_smooth(col=I("green"),method="lm",se=FALSE,na.rm=TRUE) + geom_smooth(col=I("red"),se=FALSE,na.rm=TRUE, span = 0.3) + theme_bw() + coord_cartesian(ylim=c(-0.5, 0.5))



# still a lot of noise. Let's plot a trading strategy! In particular, define 20 mutual funds that each year invest in the 
# i'th vingtile of B/M. Assume the mutual funds equal-weight the stocks each year.

# First, define the categorical variable from 1 - 20 of which B/M vingtile a stock belongs to
StockRetAcct_DT[,bm_vingtile:=cut(StockRetAcct_DT$lnBM,breaks=quantile(StockRetAcct_DT$lnBM,probs=c(0:20)/20,na.rm=TRUE), labels=FALSE)]

# get the mean excess return by vingtile
EW_BM_MutualFunds <- StockRetAcct_DT[,list(MeanExRet = mean(ExRet)), by = bm_vingtile]

# Now, plot based on the vingtile variable. This is the average excess return in the same data period. 
qplot(bm_vingtile, MeanExRet,data = EW_BM_MutualFunds, col=I("blue"), na.rm = TRUE, main = "Log Book-to-Market vs. Excess Returns") + geom_smooth(col=I("red")) + theme_bw()


# issue: This isn't really a tradeable strategy. The bins are defined on the whole sample, i.e. using forward-looking
# information. Also, there may be years where no stocks are in the extreme bins. Finally, note that we are
# implicitly weighting later years more since mean is taken across all firms and years and since there are many more
# firms in the late 1990s than early 1980s.

# We want to do the vingtile-sort each year, and we want the mutual funds to expand the set of stocks as in 
# the real world -- as they become available

# due to winsorizing of original data, we add a tiny amount of noise (jitter) to lnBM before creating by year vingtiles
# this is to avoid ties in the quantile sorts
StockRetAcct_DT[,lnBM:=jitter(lnBM, amount = 0)]

# loop through the years in the data base
for (i in 1981:2014)
{
  StockRetAcct_DT[year == i,bm_vingtile_yr:=cut(StockRetAcct_DT[year == i,]$lnBM,breaks=quantile(StockRetAcct_DT[year == i,]$lnBM,probs=c(0:20)/20,na.rm=TRUE), include.lowest=TRUE, labels=FALSE)]  
}

# now, we need to take the mean in a clever way to account for the changing number of firms
# first take mean by year and vingtile
EW_BM_MutualFunds_yr <- StockRetAcct_DT[,list(MeanExRetYr = mean(ExRet)), by = list(bm_vingtile_yr, year)]
# then average across years
EW_BM_MutualFunds_yr <- EW_BM_MutualFunds_yr[,list(MeanExRet = mean(MeanExRetYr)), by = bm_vingtile_yr]
qplot(bm_vingtile_yr, MeanExRet,data = EW_BM_MutualFunds_yr, col=I("blue"), na.rm = TRUE, main = "Book-to-Market bins vs. Excess Returns") + geom_smooth(col=I("red")) + theme_bw()


# how about splitting the data by size. Perhaps there is an interaction there?
# let's see if there is a difference between running these strategies for S&P500 firms (500 largest firms) 
# versus the remaining firms.

# first, create a rank variable based on market equity of the firms
for (i in 1981:2014)
{
  StockRetAcct_DT[year == i,size_rank := rank(-lnME)]  
}

# generate True/False dummy variable for whether in top 500 largest stocks or not
StockRetAcct_DT[,LargeStock := (size_rank<501)]  
# value
EW_BM_MutualFunds_yr_Large <- StockRetAcct_DT[LargeStock==TRUE,list(MeanExRetYr = mean(ExRet)), by = list(bm_vingtile_yr, year)]
EW_BM_MutualFunds_yr_Small <- StockRetAcct_DT[LargeStock==FALSE,list(MeanExRetYr = mean(ExRet)), by = list(bm_vingtile_yr, year)]
# then average across years
EW_BM_MutualFunds_yr_Large <- EW_BM_MutualFunds_yr_Large[,list(MeanExRet = mean(MeanExRetYr)), by = bm_vingtile_yr]
EW_BM_MutualFunds_yr_Small <- EW_BM_MutualFunds_yr_Small[,list(MeanExRet = mean(MeanExRetYr)), by = bm_vingtile_yr]
qplot(bm_vingtile_yr, MeanExRet,data = EW_BM_MutualFunds_yr_Small, col=I("blue"), ylim = c(0.025,0.125), na.rm = TRUE, main = "EW Small Firm Book-to-Market bins vs. Excess Returns") + geom_smooth(col=I("red")) + theme_bw()


# Next, let's create value-weights within the portfolio to increase tradeability. Note, we equal-weight
# across years for each portfolio, value-weight between firms within year

# value=weight
EW_BM_MutualFunds_yr_Large <- StockRetAcct_DT[LargeStock==TRUE,list(MeanExRetYr = weighted.mean(ExRet, MEwt)), by = list(bm_vingtile_yr, year)]
EW_BM_MutualFunds_yr_Small <- StockRetAcct_DT[LargeStock==FALSE,list(MeanExRetYr = weighted.mean(ExRet, MEwt)), by = list(bm_vingtile_yr, year)]
# then average across years
EW_BM_MutualFunds_yr_Large <- EW_BM_MutualFunds_yr_Large[,list(MeanExRet = mean(MeanExRetYr)), by = bm_vingtile_yr]
EW_BM_MutualFunds_yr_Small <- EW_BM_MutualFunds_yr_Small[,list(MeanExRet = mean(MeanExRetYr)), by = bm_vingtile_yr]
qplot(bm_vingtile_yr, MeanExRet,data = EW_BM_MutualFunds_yr_Small, col=I("blue"), ylim = c(0.025,0.125), na.rm = TRUE, main = "VW Small Firm Book-to-Market bins vs. Excess Returns") + geom_smooth(col=I("red")) + theme_bw()
qplot(bm_vingtile_yr, MeanExRet,data = EW_BM_MutualFunds_yr_Large, col=I("blue"), ylim = c(0.025,0.125), na.rm = TRUE, main = "VW Large Firm Book-to-Market bins vs. Excess Returns") + geom_smooth(col=I("red")) + theme_bw()





