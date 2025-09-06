# code snippets topic 8

setwd("C:/Users/llochsto/Dropbox/Data Analytics/Data")
#setwd("D:/llochsto/Dropbox/Data Analytics/DAML_2017/Data")

# Clean workspace
rm(list = ls())

# run asymptotic pca
require(MTS)
require(foreign)
require(data.table)
require(ggplot2)

# Download data and set as data.table
Factors_DT <- as.data.table(read.dta("French_Factor_Data_2012_2016.dta"))

# run asymptotic pca
French_Factors = apca(Factors_DT[, !"date"],2)
qplot(1:60,French_Factors$sdev, xlab = "factor", ylab="St.Dev. of factor",main="APCA factor standard deviations")


# return factors extracted from apca
apca_factors <- data.table(date = 1:60, factor1 = French_Factors$factors[,1], factor2 = French_Factors$factors[,2])
ggplot(data = apca_factors, aes(x = date, y = factor1)) + geom_line(color = "red") + 
  geom_line(data = apca_factors, aes(x = date, y = factor2), color = "blue") + theme_bw() +
  xlab('Months') + ylab('Factor Returns') 


# plot factor loadings
apca_loadings <- data.table(BaseAsset = 1:138, loadings1 = French_Factors$loadings[,1], loadings2 = French_Factors$loadings[,2])
ggplot(data = apca_loadings, aes(x = BaseAsset, y = loadings1)) + geom_bar(stat = "identity", alpha = 0.5, fill = "red", color = "red") + 
  geom_bar(data = apca_loadings, aes(x = BaseAsset, y = loadings2), stat = "identity", fill = "blue", alpha = 0.5, color = "blue") + theme_bw() +
  xlab('Base Assets') + ylab('Factor loadings') 




# K-means clustering

# Clean workspace
rm(list = ls())

# K-means clusters for understanding data
require(quantmod)
require(ggplot2)
Sys.setenv(TZ="GMT")
getSymbols('SPY',from='2000-01-01')

# plot density of daily SPY returns
x=remove_missing(data.frame(d=index(Cl(SPY)),return=as.numeric(Delt(Cl(SPY)))))
ggplot(x,aes(return))+stat_density(colour="steelblue", size=2, fill=NA)+xlab(label='Daily returns')

nasa=tail(cbind(Delt(Op(SPY),Hi(SPY)),Delt(Op(SPY),Lo(SPY)),Delt(Op(SPY),Cl(SPY))),-1)
colnames(nasa)=c('High','Low','Close')


# choose 5 clusters
kmeanObject=kmeans(nasa,5,iter.max=10)
kmeanObject$centers
autocorrelation=head(cbind(kmeanObject$cluster,lag(as.xts(kmeanObject$cluster),1)),-1)
autocorrelation <- remove_missing(autocorrelation)

# show autocorrelations as transition counts
xtabs(~autocorrelation[,1]+(autocorrelation[,2]))

# get autocorrelations as percentages
y=apply(xtabs(~autocorrelation[,1]+(autocorrelation[,2])),1,sum)
x=xtabs(~autocorrelation[,1]+(autocorrelation[,2]))

z=x
for(i in 1:5)
{
  z[i,]=(x[i,]/y[i])
}
round(z,2)

#optimal number of clusters
wss = (nrow(nasa)-1)*sum(apply(nasa,2,var))
for (i in 2:15) wss[i] = sum(kmeans(nasa, centers=i)$withinss)
wss=(data.frame(number=1:15,value=as.numeric(wss)))
# plot within-group sum of squares
ggplot(wss,aes(number,value))+geom_point()+
  xlab("Number of Clusters")+ylab("Within groups sum of squares")+geom_smooth()


