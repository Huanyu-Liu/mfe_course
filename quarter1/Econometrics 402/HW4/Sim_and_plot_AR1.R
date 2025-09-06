#
# simulate from AR(1)
#
beta0=0
beta1=-.8
T=100
sigma=.3
simar1(beta0,beta1,sigma,T)

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