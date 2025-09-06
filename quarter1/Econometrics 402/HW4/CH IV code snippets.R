#
# Ch IV code snippets
#
library(DataAnalytics)

data(lmich_yr)
Box.test(rnorm(100),type="Ljung",lag=20)
Box.test(lmich_yr$Level,type="Ljung",lag=20)
Box.test(lm(lmich_yr$Level~back(lmich_yr$Level))$res,type="Ljung",lag=20)

library(quantmod)
#
# get GDP
#
getSymbols("GDPC1",src = "FRED")            
#      GDP, qtrly, seasonally adjusted, 2009 $
#      note: this is an xts object not a ts object!
lnGDP=log(GDPC1)
plot(lnGDP,las=2,cex.axis=.75,main="log GDP in 2009$, SA") 
#      make labels perp to axis and reduce size
trend=1:length(lnGDP)
out=lm(lnGDP~trend)
lines(out$fitted,col="red")
acf(out$res)

acf(diff(lnGDP),na.action=na.omit)

#
# let's fit AR(1) 
#
lnGDP=as.vector(lnGDP)
nstep=50
out.ar=lm(lnGDP~back(lnGDP))
pred.ar=double(nstep+1)
pred.ar[1]=lnGDP[length(lnGDP)]
for(i in 1:nstep){pred.ar[i+1]=out.ar$coef[1]+out.ar$coef[2]*pred.ar[i]}
#
# now fit ARIMA(1,1,0)
#
out.arima=lm(diff(lnGDP)~back(diff(lnGDP)))
nstep=50
out.ar=lm(lnGDP~back(lnGDP))
pred.arima=double(nstep+1)
pred.arima[1]=lnGDP[length(lnGDP)]-lnGDP[length(lnGDP)-1]
for(i in 1:nstep){pred.arima[i+1]=out.arima$coef[1]+out.arima$coef[2]*pred.arima[i]}
pred.arima=pred.arima[-1]
pred.arima=lnGDP[length(lnGDP)]+cumsum(pred.arima)


plot(ts(c(lnGDP,pred.arima),start=c(1947,1),freq=4),type="n",ylab="lnGDP")
lines(ts(lnGDP,start=c(1947,1),freq=4),col="blue",lwd=2)
lines(ts(pred.arima,start=c(2013,4),freq=4),col="red",lty=2,lwd=2)
lines(ts(pred.ar,start=c(2013,4),freq=4),col="green",lty=3,lwd=2)
legend("topleft",legend=c("AR(1)","ARIMA(1,1,0)"),fill=c("green","red"))


#
# get Google stock quotes
#
getSymbols("GOOG")
lnGclose=log(GOOG[,4])
par(mfrow=c(2,1))
plot(lnGclose)
plot(diff(lnGclose))
par(mfrow=c(1,2))
acf(lnGclose)
acf(diff(lnGclose),na.action=na.omit)
Box.test(diff(lnGclose),type="Ljung",lag=20)

#
# simulate from AR(1)
#
beta0=0
beta1=.8
T=100
sigma=.3
y=simar1(beta0,beta1,sigma,T)

simar1=function(beta0,beta1,sigma,T){
  mu=beta0/(1-beta1)
  y=double(T)
  y[1]=mu
  for(t in 2:T){y[t]=beta0+beta1*y[t-1]+rnorm(1,sd=sigma)}
  plot(y,type="n",ylab="",xlab="time")
  #
  # color in positive and negative parts
  #
  xint=function(x,y){b=(y[2]-y[1])/(x[2]-x[1]);a=y[1]-b*x[1];xint=(mu-a)/b}
  yoldamu=y[1]>mu
  for(t in 2:T){
    yamu=y[t]>mu
    if(yoldamu)         # here y_t01 is above mean
    {if(yamu)          # here y_t above mean
    {polygon(x=c(t-1,t-1,t,t),y=c(mu,y[t-1],y[t],mu),col="green",lty=0,border="white")}
     else              # here y_t below mean
     {xintercept=xint(c(t-1,t),c(y[t-1],y[t]))
      polygon(x=c(t-1,t-1,xintercept),y=c(mu,y[t-1],mu),col="green",lty=0,border="white")
      polygon(x=c(xintercept,t,t),y=c(mu,mu,y[t]),col="red",lty=0,border="white")}}
    else                 # here y_t-1 is below the mean
    {if(yamu)           # here y_t is above the mean
    {xintercept=xint(c(t-1,t),c(y[t-1],y[t]))
     polygon(x=c(t-1,xintercept,t-1),c(mu,mu,y[t-1]),col="red",lty=0,border="white")
     polygon(x=c(xintercept,t,t),y=c(mu,y[t],mu),col="green",lty=0,border="white")}
     else               # here y_t is below the mean
     {polygon(x=c(t-1,t,t,t-1),y=c(mu,mu,y[t],y[t-1]),col="red",lty=0,border="white")}
    }
    yoldamu=yamu
  }
  points(y,pch=20,col="blue")
  lines(y)
  abline(h=mu)
  title(paste("AR(1), beta_0 =",beta0,", beta_1 =",beta1,sep=""))
  invisible(y)
}


# retrieve Google Trends data
#
data(gtAuto)
#from http://www.census.gov/retail/marts/www/timeseries.html
gtAuto$y=log(gtAuto$sales) #log of sales (sales are in millions of dollars)

qplot(Period,y,data=gtAuto,main="log Auto Sales",xlab="",ylab="") + 
  geom_point(col=I("red"),size=I(3)) + 
  geom_line(col=I("blue"),linetype="dashed")+theme_bw()

acf(gtAuto$y)
lmout_base=lm(y~back(y)+back(y,12),data=gtAuto)
lmSumm(lmout_base)
acf(lmout_base$resid)
#
# now add in Google Trends indices of search activity for terms related to
# auto sales
#
lmout_trends=lm(y~back(y)+back(y,12)+suvs+insurance,data=gtAuto)
lmSumm(lmout_trends)
#
# perform one-step ahead prediction
gtAuto$y.lag1 = back(gtAuto$y) #1 period lagged sales
gtAuto$y.lag12 = back(gtAuto$y,12) #12 period lagged sales
n = length(gtAuto$y)
k=17 #start with only the first 17 observations and then add one by one
gtAuto$Actual=gtAuto$Baseline=gtAuto$Trends=rep(NA,n) #initialize columns

for (t in k:(n-1)) {
  # roll forward the regressions to predict one-step ahead
  reg1 = lm(y~y.lag1+y.lag12,data=gtAuto[1:t,])        
  reg2 = lm(y~y.lag1+y.lag12+suvs+insurance, data=gtAuto[1:t,]) 
  t1 = t+1
  gtAuto$Actual[t1]=gtAuto$y[t1]
  gtAuto$Baseline[t1]=predict(reg1,newdata=gtAuto[t1,])
  gtAuto$Trends[t1]=predict(reg2,newdata=gtAuto[t1,])
}

# remove rows and cols in dataframe for plotting
z = gtAuto[(k+1):n,c("Period","Actual","Baseline","Trends")]

qplot(Period,y=Actual, color = "Actual", data=z,geom="line",
      main="Motor Vehicles and Parts",ylab="log(sales)") +
  geom_line(aes(y = Baseline, color = "Baseline")) +
  geom_line(aes(y = Trends, color = "Trends")) +
  theme_bw()  

mae1 = mean(abs(exp(z$Actual)-exp(z$Baseline))/exp(z$Actual))
mae2 = mean(abs(exp(z$Actual)-exp(z$Trends))/exp(z$Actual))
mae2/mae1-1

#
# compare to insample fit
#
ActualSales=gtAuto$sales[13:length(gtAuto$sales)]
mae1_insam=mean(abs(ActualSales-exp(lmout_base$fitted)))
mae2_insam=mean(abs(ActualSales-exp(lmout_trends$fitted)))
mae2_insam/mae1_insam-1
