library("xts")
mydata = read.csv("feds200628.csv")
mydata$X = as.Date(mydata$X)
rownames(mydata) = mydata[,1]
mydata$X = NULL
start_yield = mydata["1983-12-30",c(1:30)]
#plot(x=c(1:30),start_yield,pch=20,type="l",lty=1)
end_yield = mydata["2018-06-29",c(1:30)]
#plot(x=c(1:30),y=end_yield,col="red",lty=2,add=TRUE)
start_yield
start_y2 = start_yield[[2]] / 100
start_y10 = start_yield[[10]] / 100
modified_D2 = 2 / (1 + start_y2)
modified_D10 = 10 / (1 + start_y10)
start_p2 = 100 / (1 + start_y2)^2
start_p10 = 100 / (1 + start_y10)^10

A = matrix(data=c((modified_D2 * start_p2), (modified_D10 * start_p10), -start_p2, start_p10), nrow=2, ncol=2, byrow=TRUE)    
b = matrix(data=c(0,1000000 /0.1),nrow=2, ncol=1, byrow=FALSE)
result = solve(A,b)
result
ini_cash = abs(result[1] * start_p2) - abs(result[2] * start_p10) + 1000000
mydata = as.xts(mydata)
mydata = mydata["1983-12-30/2018-06-29"]
NSS_model = function(parameters,t){
  beta0 = parameters[1]
  beta1 = parameters[2]
  beta2 = parameters[3]
  beta3 = parameters[4]
  tau1 = parameters[5]
  tau2 = parameters[6]
  t1 = t/tau1
  t2 = t/tau2
  e1 = exp(-t1)
  e2 = exp(-t2)
  return(beta0 + beta1 * (1 - e1)/t1 + beta2 * ((1 - e1)/t1 - e1) + beta3 * ((1 - e2) / t2 - e2))
}
weekly = endpoints(mydata,"weeks", k=1)
weekly[3]

firstweek = as.vector(mydata[weekly[3],])
parameter = firstweek[(length(firstweek) - 5):length(firstweek)]
t10 = 9+358/365
t2 = 1 + 358/365
y10 = NSS_model(parameters= parameter,t=t10)/100
price10 = 100 / (1 + y10)^(t10)
price10
y2 = NSS_model(parameters = parameters, t = t2)/100
price2 = 100/(1+y2)^(t2)
week_y = ((1 + start_y2) ^ 2 / (1 + y2) ^ t2) ^ (365/7) - 1

startweek = as.vector(mydata[weekly[2],])
parameter = startweek[(length(startweek) - 5):length(startweek)]
t_week = 7 / 365
y_w = NSS_model(parameters= parameter,t=t_week)/100
interest = ini_cash * (1 + y_w/(365/7))^(t_week) - ini_cash


cash = ini_cash + interest + result[2,1] * price10 + result[1,1]* price2

week_return = cash /1000000 - 1

weekly = weekly[-1]
initial_capital = 1000000
myfunction = function(initial_capital){
  start_yield = mydata["1983-12-30",c(1:30)]
  end_yield = mydata["2018-06-29",c(1:30)]
  for (i in weekly){
    startweek = as.vector(mydata[i,])
    sw_parameter = startweek[(length(startweek) - 5):length(startweek)]
    t_week = 7 / 365
    y_w = NSS_model(parameters= parameter,t=t_week)/100
    
  }
}


  