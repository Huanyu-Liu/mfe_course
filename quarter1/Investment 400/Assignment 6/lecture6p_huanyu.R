mysoln = list(student = c("Huanyu Liu", "Hyeuk Jung", "Jiaqi Li", "Xichen Luo"))
library("quantmod")
library("xts")
getSymbols("MSFT",from="1989-12-29",to="2018-09-29")
getSymbols("INTC",from="1989-12-29",to="2018-09-29")
MSFT = as.xts(MSFT)
INTC = as.xts(INTC)
result = list()

myfunc = function(data){
  week_end = endpoints(data,"weeks",k=1)[-1]
  weekly_return = vector()
  for (i in 1:(length(week_end) - 1)){
    weekly_return[i] = data[[week_end[i+1],6]] / data[[week_end[i],6]] - 1
  }
  mean = mean(weekly_return)
  sd = sd(weekly_return)
  annualized_mean = (1 + mean)^52 -1
  annualized_sd = sd * sqrt(52)
  result[[1]] = weekly_return
  result[[2]] = c(mean,sd,annualized_mean,annualized_sd)
  return(result)
}
stats_ms = myfunc(MSFT)
stats_it = myfunc(INTC)
mysoln[["Q1"]] = list("intel" = stats_it[[2]],"microsoft" = stats_ms[[2]])

R_f = 0.01
A = 4
w_it = (stats_it[[2]][3] - R_f)/(A * stats_it[[2]][4]^2)
w_f_it = 1 - w_it

w_ms = (stats_ms[[2]][3] - R_f)/(A * stats_ms[[2]][4]^2)
w_f_ms = 1 - w_ms

mysoln[["Q2"]] = list("intel-RF" = c(w_it, w_f_it),"microsoft-RF" = c(w_ms,w_f_ms))


r_it = w_it * stats_it[[2]][3] + w_f_it * R_f
r_ms = w_ms * stats_ms[[2]][3] + w_f_ms * R_f
r_it - r_ms

#choose microsoft, because of the higher return of the portforlio
mysoln[["Q3"]] = c("microsoft",r_it - r_ms)
correlation = cor(stats_ms[[1]],stats_it[[1]])
return_func = function(weight){
  return(weight * stats_it[[2]][3] + (1 - weight) * stats_ms[[2]][3])
}
sd_func = function(weight){
  return(sqrt(weight^2 * stats_it[[2]][4]^2 + (1 - weight)^2 * stats_ms[[2]][4]^2 + 2 * weight * (1 - weight) * correlation * stats_it[[2]][4] * stats_ms[[2]][4]))
}

returns = return_func(c(1:10000)/10000)
sds = sd_func(c(1:10000)/10000)
plot(x=sds,y=returns,type = "l")
which.min(sds)
min_var_it_weight = 3267 / 10000
min_var_ms_weight = 1 - min_var_it_weight
min_var_return = return_func(min_var_it_weight)
mysoln[["Q4"]] = list("intel_weight" = min_var_it_weight, "microsoft_weight" = min_var_ms_weight, "min_return" = min_var_return, "min_var" = min(sds))
mysoln