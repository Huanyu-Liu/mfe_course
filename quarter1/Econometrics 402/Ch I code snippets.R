#
# M 237Q
# Chapter I code snippets  
#
#
# install course package and ggplot2
install.packages("devtools")
install.packages("ggplot2")
library(devtools)
library(ggplot2)


install_bitbucket("perossichi/DataAnalytics")


library(DataAnalytics)

#
# reshape Vanguard data using reshape2
#
# or use the reshape package
#
library(reshape2)
data(Vanguard)
Van=Vanguard[,c(1,2,5)]   # grab relevant cols
V_reshaped=dcast(Van,date~ticker,value.var="mret")

#
# now let's plot mean std deviation.
#
mat=descStat(V_reshaped)
#
SD=mat[,3]
MEAN=mat[,1]
funds=rownames(mat)
plot(SD,MEAN,pch=20,col=rainbow(length(funds)),cex=2)
text(SD,MEAN,funds,col=rainbow(length(funds)),pos=1)
#
# do it in ggplot2
#
qplot(SD,MEAN,geom=c("point","text"),label=funds,vjust=-1,xlab="StdDev",ylab="Mean")
#
# or color
#
qplot(SD,MEAN,xlab="StdDev",ylab="Mean") +
  geom_point(col=rainbow(length(funds)),size=4) +
  geom_text(label=funds,col=rainbow(length(funds)),vjust=-1)

#
# let's merge in market indices and plot
#
data(marketRf)
Van_mkt=merge(V_reshaped,marketRf,by="date")
with(Van_mkt,
     plot(vwretd,VGHCX,pch=20,col="blue")
)
     
out=lm(VGHCX~vwretd,data=Van_mkt)
abline(out$coef,col="red",lwd=2)
with(Van_mkt,
     points(mean(vwretd),mean(VGHCX),pch=20,cex=2,col="green")
)
#
# compute intercept for line with slope b going thru (mean(x),mean(y))
#
badslope=-.5
int=mean(Van_mkt$VGHCX)-badslope*mean(Van_mkt$vwretd)
abline(c(int,badslope),lty=2,lwd=2,col="green")

#
# plot with marginal distribution
#
with(Van_mkt,
       plot(vwretd,VGHCX,pch=20,col="blue")
)
rug(Van_mkt$VGHCX,side=2,col="red")
qvec=quantile(Van_mkt$vwretd,probs=c(.2,.4,.6,.8))
abline(v=qvec,lty=2,lwd=2,col="green")
#
# plot each of conditionals
#
# slice using "cut" function
cat=cut(Van_mkt$VGHCX,breaks=c(min(Van_mkt$VGHCX),qvec,max(Van_mkt$VGHCX)))
cat=c(as.numeric(cat),rep(6,length(cat)))
cat=factor(cat,
           labels=c("0-20","21-40","41-60","61-80","81-100","all"))
# repeat whole data on bottom to get marginal
VGHCXbig=rep(Van_mkt$VGHCX,2)
boxplot(VGHCXbig~cat,col="green")

#
# let's simulate sampling distribution of b_1
#
simreg=function(beta0,beta1,sigma,x){
  y=beta0+beta1*x+rnorm(length(x),sd=sigma)
}
#

x=Van_mkt$vwretd
beta0=.006
beta1=.75
sigma=0.025
nsample=10000
bsim=double(nsample)
for(i in 1:nsample){
  y=simreg(beta0,beta1,sigma,x)
  bsim[i]=lm(y~x)$coef[1]
}
hist(bsim,breaks=40,col="magenta")

#
# example of hypothesis testing
#
t=(.743-1)/.030130
pvalue=2*pt(-abs(t),df=347)
pvalue

#
# predict new values from HCX on VWRET regression
#
library(reshape2)
data(Vanguard)
Van=Vanguard[,c(1,2,5)]   # grab relevant cols
V_reshaped=dcast(Van,date~ticker,value.var="mret")
data(marketRf)
Van_mkt=merge(V_reshaped,marketRf,by="date")

out=lm(VGHCX~vwretd,data=Van_mkt)
predict(out,new=data.frame(vwretd=.10),int="prediction")
