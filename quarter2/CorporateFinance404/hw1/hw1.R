library("xts")
library("quantmod")
df = read.csv("sp500.csv")
getSymbols("DGS3MO",from = as.Date('1970-01-01'), to = as.Date('2015-05-31'),src = "FRED")
getSymbols("DGS1", src="FRED")
getSymbols("DGS5", src="FRED")
rownames(df) = as.Date(as.character(df$caldt),"%Y%m%d")
df$caldt = NULL
sp500 = as.xts(df)
daily_returns = df$vwretd

geometric_mean = function(v){
  geo = (prod(v + 1)^(1/length(v)) - 1)
  return(geo)
}

year_end = endpoints(sp500,"years",k=1)
month_end = endpoints(sp500,"months",k =1)

arith_daily_sp500 =  mean(df$vwretd)
geo_daily_sp500 = geometric_mean(daily_returns)

monthly_return = vector()
annual_return = vector()
year_5_return = vector()

for (i in c(1:(length(month_end) - 1))){
  monthly_return[i] = prod(daily_returns[(month_end[i] + 1) : month_end[i+1]] + 1) - 1
}

arith_monthly_sp500 = mean(monthly_return)
geo_monthly_sp500 = geometric_mean(monthly_return)

for (i in c(1:(length(year_end) - 1))){
  annual_return[i] = prod(daily_returns[(year_end[i] + 1):year_end[i+1]] + 1) - 1
}

arith_annual_sp500 = mean(annual_return)
geo_annual_sp500 = geometric_mean(annual_return)

for (i in c(1:(length(year_end) %/% 5))){
  year_5_return[i] = prod(daily_returns[(year_end[i*5 - 4] + 1):(year_end[i*5 + 1])] + 1) - 1
}
arith_5year_sp500 = mean(year_5_return)
geo_5year_sp500 = geometric_mean(year_5_return)

