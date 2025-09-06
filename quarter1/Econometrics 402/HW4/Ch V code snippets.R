#
# ch V code snippets
#
library(DataAnalytics)

#
# let's check for hidden extrapolation
#
data(multi)
fullout=lm(Sales~p1+p2,data=multi)
X=cbind(c(rep(1,nrow(multi))),multi$p1,multi$p2)
H=X%*%solve(crossprod(X))%*%t(X)
x_f=c(1,8,6)
x_f=as.matrix(x_f,ncol=1)
hat_f=t(x_f)%*%chol2inv(chol(crossprod(X)))%*%x_f
hat_f
hist(hatvalues(fullout),main="",col="magenta")

#
# check residual diagnostics
#
library(reshape2)
data(Vanguard)
Van=Vanguard[,c(1,2,5)]   # grab relevant cols
V_reshaped=dcast(Van,date~ticker,value.var="mret")
data(marketRf)
Van_mkt=merge(V_reshaped,marketRf,by="date")
residualPlots(lm(VSMAX~vwretd,data=Van_mkt))
#
# check hedge fund return performance
#
data(HFRI)
HFRImkt<-merge(HFRI,marketRf,by="date")
lmSumm(lm(EHQD~vwretd+back(vwretd),data=HFRImkt),L=10,HAC=TRUE)
acf(lm(EHQD~vwretd,data=HFRImkt)$res)

#
# let's bootstrap a regression
#
library(reshape2)
data(Vanguard)
Van=Vanguard[,c(1,2,5)]   # grab relevant cols
V_reshaped=dcast(Van,date~ticker,value.var="mret")
data(marketRf)
Van_mkt=merge(V_reshaped,marketRf,by="date")


Reg_data=data.frame(Van_mkt$vwretd,Van_mkt$VGENX)
Reg_data=na.omit(Reg_data)
colnames(Reg_data)=c("GENX","VW")
N=nrow(Reg_data)
regression=summary(lm(GENX~VW,data=Reg_data))
lsq.slope=regression$coef[2,1]
lsq.slope.stderr=regression$coef[2,2]
#
# add in exp(slope)
#

B=10000
BS_coefs=matrix(0,nrow=B,ncol=2)
for(b in 1:B){
  BS_sample=Reg_data[sample(1:N,size=N,replace=TRUE),]
  sum=summary(lm(GENX~VW,data=BS_sample))
  BS_coefs[b,]=sum$coef[,1]
  }

par(mfrow=c(1,2))
hist(BS_coefs[,2],breaks=40,col="magenta",xlab="",ylab="",
     main="Bootstrap Dist of Slope",xlim=c(.3,.60),
     sub=paste("bootstrap std error = ",round(sd(BS_coefs[,2]),digits=3)))
summary(BS_coefs[,2])
hist(rnorm(B,mean=lsq.slope,sd=lsq.slope.stderr),breaks=40,col="magenta",
     xlab="",ylab="",main="Normal Dist of Slope",xlim=c(.3,.60),
     sub=paste("Least Squares Std Error = ",round(lsq.slope.stderr,digits=3)))
#
# let's construct a Bootstrap CI for the slope coef
#
int=quantile(BS_coefs[,2],probs=c(.975,.025))
CI.pivotal.bootstrap = c(2*lsq.slope-int[1],2*lsq.slope-int[2])
CI.normal.theory=c(lsq.slope+qt(.025,df=347)*lsq.slope.stderr,
                   lsq.slope+qt(.975,df=347)*lsq.slope.stderr)
CI.mat=rbind(CI.normal.theory,CI.pivotal.bootstrap)
colnames(CI.mat)=c("Lower","Upper")
rownames(CI.mat)=c("Normal Theory","Bootstrap Pivotal")
print(CI.mat)

#
# bootstrap exp(slope)
#
par(mfrow=c(1,1))
hist(exp(BS_coefs[,2]),breaks=40,col="magenta",xlab="",ylab="",
     main="Bootstrap Distribution of exp(Slope)",
     sub=paste("bootstrap std error = ",round(sd(exp(BS_coefs[,2])),digits=3)))
abline(v=exp(lsq.slope),lwd=3,col="yellow")
#
# and pivot CI
int=quantile(exp(BS_coefs[,2]),probs=c(.975,.025))
CI.pivotal.bootstrap = c(2*exp(lsq.slope)-int[1],2*exp(lsq.slope)-int[2])
names(CI.pivotal.bootstrap)=c("Lower","Upper")
CI.pivotal.bootstrap
  

