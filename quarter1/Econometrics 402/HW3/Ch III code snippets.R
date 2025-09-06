#
# Chapter III code snippets
# let's play around with linear algebra
#
X=matrix(rnorm(30000),ncol=3)  # fill a matrix with numbers, matrix is 10,000 x 3
dim(X)
dim(t(X))
crossprod(X)
chol(crossprod(X))
#
#
library(reshape2)
data(Vanguard)
Van=Vanguard[,c(1,2,5)]   # grab relevant cols
V_reshaped=dcast(Van,date~ticker,value.var="mret")
data(riskFactors)
Van_risk=merge(V_reshaped,riskFactors,by="date")
#
# let's do it using matrix algebra
y=Van_risk$VWNFX
X=cbind(rep(1,length(y)),Van_risk$RmRf,Van_risk$SMB,Van_risk$HML)
# remove missing from X and y
X=X[-which(is.na(y)),]
y=y[-which(is.na(y))]
b=chol2inv(chol(crossprod(X)))%*%crossprod(X,y)
#
# compute std errors
#
e=y-X%*%b
ssq=sum(e*e)/(length(y)-ncol(X))
Var_b=ssq*chol2inv(chol(crossprod(X)))
std_err=sqrt(diag(Var_b))
names(std_err)=c("intercept","RmRf","SMB","HML")

#
# simulate from MVN
#

k=2
n=1000
mu=c(rep(1,k))
Sigma=matrix(c(1,.8,.8,1),ncol=2)
U=chol(Sigma)
Z=matrix(rnorm(n*k),ncol=n)
Y=crossprod(U,Z)+mu  # Y is k x n
Y=t(Y)    #Y is n x k now
# or
dim(Z)=c(n,k)
Y=Z%*%U + mu

#
#
# chi-squared distribution
#
x=seq(from=.001,to=12,len=100)
plot(x,dchisq(x,df=3),lwd=2,col="blue",type="l",xlab="",ylab="",
     main="Chi-Squared with 3 df")



#
# let's do a market timing test on Vanguard data
#
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
data(marketRf)
Van_mkt=merge(V_reshaped,marketRf,by="date")
mkt_up=ifelse(Van_mkt$vwretd>0,1,0)
Van_mkt$upvw=mkt_up*Van_mkt$vwretd
Van_mkt$dwnvw=(1-mkt_up)*Van_mkt$vwretd
mkt_timing=lm(VWNFX~upvw+dwnvw,data=Van_mkt)
lmSumm(mkt_timing)
#
# now compute GLH F stat
#
R=matrix(c(0,1,-1),byrow=TRUE,nrow=1)
r=c(0)
X=cbind(c(rep(1,nrow(Van_mkt))),Van_mkt$upvw,Van_mkt$dwnvw)
b=as.vector(mkt_timing$coef)
QFmat=chol2inv(chol(crossprod(X)))
QFmat=R%*%QFmat%*%t(R)
Violation=R%*%b-matrix(r,ncol=1)
fnum=t(Violation)%*%chol2inv(chol(QFmat))%*%Violation
n_minus_k = nrow(Van_mkt)-length(b)
fdenom=nrow(R)*sum(mkt_timing$resid**2)/n_minus_k
f=fnum/fdenom
f
pvalue=1-pf(f,df1=nrow(R),df2=n_minus_k)
pvalue

