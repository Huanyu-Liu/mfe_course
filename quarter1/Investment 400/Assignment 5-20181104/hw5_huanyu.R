library("xts")
sp500 = read.csv("p5-sp500.csv")
sp500$caldt = as.Date(as.character(sp500$caldt),"%Y%m%d")
rownames(sp500) = sp500[,1]
sp500$caldt = NULL
sp500 = as.xts(sp500)

year_end = endpoints(sp500,"years",k=1)
month_end = endpoints(sp500,"months",k =1)

annual_return = vector()
month_return = vector()
year5_return = vector()

daily_returns = sp500["1972-01-03/2017-12-30"]
daily_return_a = mean(daily_returns) * 365
daily_return_g = (prod(daily_returns + 1)^(1/length(daily_returns)) - 1) * 365


for (i in 1:(length(month_end) - 1 )){
  returns = as.vector(sp500[(month_end[i] + 1):month_end[i+1]])
  month_return[i] = prod(1 + returns) - 1
}

month_return_a = mean(month_return) * 12
month_return_g = (prod(month_return + 1)^(1/length(month_return)) - 1) * 12

for (i in 1:(length(year_end) - 1)){
  returns = as.vector(sp500[(year_end[i]+1):year_end[i+1]])
  annual_return[i] = prod(1 + returns) - 1
}


annual_return_a = mean(annual_return)
annual_return_g = prod(annual_return + 1)^(1/(2017 - 1972 + 1)) -1

for (i in seq(1,length(year_end) - 6, 5)){
  print(i)
  print(year_end[i])
  print(sp500[year_end[i + 5]])
  returns = as.vector(sp500[(year_end[i] + 1):year_end[i+5]])
  year5_return[i %/% 5 + 1] = prod(returns + 1) - 1
}

year5_return_a = mean(year5_return) / 5
year5_return_g = prod(1 + year5_return)^(1/length(year5_return)) / 5


daily_return_a
daily_return_g
month_return_a
month_return_g
annual_return_a
annual_return_g
year5_return_a
year5_return_g