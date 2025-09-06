#
# code snippets for Topic 3
#
# set your working directory, where you have downloaded the Stata data file
# change the below to match your folder
setwd("D:/llochsto/Dropbox/Data Analytics/Data")
# setwd("C:/Users/Lars/Dropbox (Personal)/Data Analytics/Data")

require(data.table)
require(ggplot2)

#
# an example of a logistic regression
#
# simulate data
#
set.seed(666)
FICO=rnorm(2000,mean=650,sd=200)
FICO=FICO[FICO< 851]
FICO=signif(FICO,2)
beta1=-.005
beta0=.5
pr=exp(beta0+beta1*FICO)
pr=pr/(1+pr)
Default=as.numeric(runif(length(FICO))<pr)
default=data.frame(FICO=FICO,Default=Default)

# plot data
qplot(FICO,Default,col=I("blue"),alpha=I(1/3),size=I(5),data=default)+
  theme(axis.title=element_text(size=rel(1.5)),
        axis.text=element_text(size=rel(1.25),colour=I("red")))
qplot(factor(Default),FICO,geom="boxplot",fill=I("green"),xlab="Default") +
  theme(axis.title=element_text(size=rel(1.5)),
        axis.text=element_text(size=rel(1.25),colour=I("red")))

# fit logit model
out=glm(Default~FICO,family="binomial",data=default)
summary(out)

#
# graph likelihood for coin toss
#
N=10
nhead=3
like=function(theta,nhead,N){(theta**nhead)*(1-theta)**(N-nhead)}
thetas=seq(from=0,to=1,length=50)
likeval=like(thetas,nhead,N)
qplot(thetas,likeval,geom="line",col=I("red"),ylab="likelihood",xlab="Theta") +
  geom_vline(xintercept=.3,colour=I("green")) +
  theme(axis.title=element_text(size=rel(1.5)),
        axis.text=element_text(size=rel(1.25),colour=I("blue")))

#
# maximize likelihood
#
y=default$Default
X=default$FICO
X=cbind(c(rep(1,length(y))),FICO)
loglike=function(beta,y,X){
  ind=X%*%beta
  pr=exp(ind)
  pr=pr/(1+pr)
  sum(y*log(pr)+(1-y)*log(1-pr))
}
beta = c(rep(0, 2))
mle = optim(beta, loglike, X = X, y = y, method = "BFGS", hessian = TRUE, 
            control = list(fnscale = -1))
mle$par
sqrt(diag(solve(-mle$hessian)))
summary(glm(Default~FICO,data=default,family="binomial"))


#
# graph likelihood for logit

ngrid=50
beta0=seq(from=0,to=1,length=ngrid)
beta1=seq(from=-.01,to=.0,length=ngrid)

ll=matrix(nrow=ngrid,ncol=ngrid)
for (i in 1:ngrid){
  for (j in 1:ngrid){
    beta=c(beta0[i],beta1[j])
    ll[i,j]=loglike(beta,y=y,X=X)
  }
}
image(beta0,beta1,ll,col=terrain.colors(60),zlim=c(-500,-400))
#
# compute the change in probability of defaulty given
# FICO goes from 500 to 800
#
predout=predict(out,type="response",new=data.frame(FICO=c(500,800)))
(predout[1]-predout[2])
#
# now let's to the analogue of the overall F test for the significance of the regression --
#  this is now a Chi-squared test with degrees of freedom equal to the number of regressors in the 
#  Logisic model
#
test_stat=out$null.deviance-out$deviance  # Chi-sq with k degrees of freedom
test_stat
k=out$df.null-out$df.residual
k
pvalue_chisq=1-pchisq(test_stat,df=k)
pvalue_chisq

anova(out)

#
# now let's plot data and fitted probs
#
plot(FICO,Default,type="n")
points(FICO[Default==1],Default[Default==1],pch=20,col="red")
points(FICO[Default==0],Default[Default==0],pch=20,col="green")
points(FICO,out$fitted,pch=20,col="blue")
legend("topleft",c("Default","No-Default","Fitted Prob"),fill=c("red","green","blue"))

#
# let's look at the loan data
#
require(DataAnalytics)
data(loans)
loans$default=loans$Good.Loan-1
loans$Good.Loan=NULL # remove old variable
names(loans)=c("StatChkA","Duration","CrdHist","Purpose","CrdAmt","Sav_Bnd",
               "Emply","InstallRate","Pstatus","OthrDebt",
               "Resid","Proprty","Age","OthrInstall","Housing","Ncredits","job","Nsupport",
               "Telephone","Foreign","default")
outloans_full=glm(default~.,data=loans,family="binomial")
summary(outloans_full)
# let's stripdown the names so that we can view the output on one "page"

#
# three types of variables
# Credit History
# loan characteristics
# demographics
#
# let's refit with only the factors that have significant coefficients 
#
outloans1=glm(default~StatChkA + CrdHist + Sav_Bnd + OthrDebt + OthrInstall +
                Purpose + Duration + CrdAmt + InstallRate  +
                Pstatus + Foreign,
              data=loans,family="binomial")
#
# compute analogue of partial F test
#
delta_df=summary(outloans_full)$df[1]-summary(outloans1)$df[1]
delta_df
delta_dev=outloans1$dev-outloans_full$dev
delta_dev
1-pchisq(delta_dev,df=delta_df)

#
# lets look at the fitted probabilities 
#
phat=predict(outloans1,type="response")
qplot(factor(loans$default),phat,geom="boxplot",fill=I("green"),xlab="Default") +
  theme(axis.title=element_text(size=rel(1.5)),
        axis.text=element_text(size=rel(1.25),colour=I("red")))


#
# compute a "lift" table
#
phat=predict(outloans1,type="response")
deciles=cut(phat,breaks=quantile(phat,probs=c(seq(from=0,to=1,by=.1))),include.lowest=TRUE)
deciles=as.numeric(deciles)
df=data.frame(deciles=deciles,phat=phat,default=loans$default)
lift=aggregate(df,by=list(deciles),FUN="mean",data=df) # find mean default for each decile
lift=lift[,c(2,4)]
lift[,3]=lift[,2]/mean(loans$default)
names(lift)=c("decile","Mean Response","Lift Factor")
lift


# compute a ROC curve
# define a function that creates the true and false positive rates
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

glm_simple_roc <- simple_roc(loans$default=="1", phat)
TPR <- glm_simple_roc$TPR
FPR <- glm_simple_roc$FPR

# plot the corresponding ROC curve
q <- qplot(FPR,TPR,xlab="FPR",ylab="TPR",col=I("blue"),main="ROC Curve for Logistic Regression Default Model",size=I(0.75))
# add straight 45 degree line from 0 to 1
q + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size=I(1.0)) + theme_bw()

# let's see how the unrestricted model fares
phat_full=predict(outloans_full,type="response")
glm_simple_roc_full <- simple_roc(loans$default=="1", phat_full)

glm_simple_roc <- cbind(glm_simple_roc,Model = "Restricted")
glm_simple_roc_full <- cbind(glm_simple_roc_full, Model = "Full")

New_ROC <- rbind(glm_simple_roc, glm_simple_roc_full)

q <- qplot(FPR,TPR,data = New_ROC, colour = Model, xlab="FPR",ylab="TPR",main="ROC Curve for Full and Restricted Logistic Regression Models",size=I(0.75))
q + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size=I(1.0), col = I("black")) + theme_bw()
# not much of a difference, as expected given model selection threw out insignificant variables


#
# Analyze Google data
#
col.alpha = function(color,alpha)
{
  code = col2rgb(color)/256
  rgb(code[1],code[2],code[3],alpha)
}
#

# download the google.csv dataset from CCLE
nexus_df <- read.csv("google.csv", header=TRUE)

#
# plot distribution of number of videos watched in electronics by watched.promo
qplot(num.videos.watched.in.electronics,geom="density",
      colour=factor(watched.promo),size=I(1.5), group=watched.promo,data=nexus_df) +
  theme_bw(base_size=14)

# fit treatment effect with no controls
#
fit.naive=glm(bought.nexus~watched.promo,data=nexus_df,family="binomial")
#
# naive estimate 
#
effect_naive = predict(fit.naive,new=data.frame(watched.promo=1),type="response") - 
  predict(fit.naive,new=data.frame(watched.promo=0),type="response")

#
# fit propensity score model
#
prop.fit <- nexus_df[, setdiff(names(nexus_df), 'bought.nexus')]
prop.out <- glm(watched.promo ~ ., data=prop.fit, family=binomial(logit))
summary(prop.out)
pscore=prop.out$fitted
nexus_df$pscore=pscore
#
# plot distributions of pscores
#
hist(nexus_df$pscore[nexus_df$watched.promo==1],breaks=20,col=col.alpha("red",.5),freq=FALSE,
     ylim=c(0,5),xlab="Propensity Score",ylab="",main="")
hist(nexus_df$pscore[nexus_df$watched.promo==0],breaks=20,col=col.alpha("blue",.5),
     freq=FALSE, ylim=c(0,5),add=TRUE)

#
# now add propensity score to the Naive model
#
fit.pscore = glm(bought.nexus~watched.promo+pscore,data=nexus_df,family="binomial")
summary(fit.pscore)
#
# compute effect for treated.  
#
ps_trt=nexus_df$pscore[nexus_df$watched.promo==1]
n=length(ps_trt)
effect_treated = 
  predict(fit.pscore,new=data.frame(watched.promo=c(rep(1,n)),pscore=ps_trt),
          type="response") - 
  predict(fit.pscore,new=data.frame(watched.promo=c(rep(0,n)),pscore=ps_trt),
          type="response")
mean(effect_treated)

#
# let's compute a propensity score based on all variables except electronics videos and age
#
prop.out <- glm(watched.promo ~ .- bought.nexus - num.videos.watched.in.electronics - pscore 
                - age,
                data=nexus_df, family=binomial(logit))
pscore=prop.out$fitted
nexus_df$pscore=pscore
fit.pscore = glm(bought.nexus~watched.promo+pscore,data=nexus_df,family="binomial")

#
# compute effect for treated.  
#
ps_trt=nexus_df$pscore[nexus_df$watched.promo==1]
n=length(ps_trt)
effect_treated = 
  predict(fit.pscore,new=data.frame(watched.promo=c(rep(1,n)),pscore=ps_trt),
          type="response") - 
  predict(fit.pscore,new=data.frame(watched.promo=c(rep(0,n)),pscore=ps_trt),
          type="response")
mean(effect_treated)






