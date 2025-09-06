library("readxl")
library("sandwich")
# Problem 1

# 3

plot_func = function(returns){
  plot(returns, pch=20, col="blue", ylim = c(-.3, .3), type = 'l')
  abline(h = mean(returns) + sd(returns)*4, col = "red") # r + sigma*4
  abline(h = mean(returns) - sd(returns)*4, col = "red")
}

r_simple = rnorm(600,mean = 0.008,sd = 0.063)
plot_func(r_simple)

mu = 0.012
sigma = 0.05
p = 0.15
mu_j = -0.03
sigma_j = 0.1
Bernoulli = rbinom(600,1,p)
J = Bernoulli * rnorm(600,mean = mu_j,sd = sigma_j)
r_jump = J + rnorm(600,mean = mu, sd = sigma)



skewness = function(return){
  mean = mean(return)
  sigma = sd(return)
  return(mean((return - mean)^3/sigma^3))
}

kurtosis = function(return){
  mean = mean(return)
  sigma = sd(return)
  return(mean((return - mean)^4/sigma^4))
}

plot_func(r_jump)

mean_jump = mean(r_jump)
sd_jump = sd(r_jump)
sk_jump = skewness(r_jump)
kurt_jump = kurtosis(r_jump)




# Problem 2

file <- "C:/Users/LIUH0/Documents/mfe courses/quarter2/EmpiricalMethod406/hw1/DBV.xlsx"
gspc_file <- "C:/Users/LIUH0/Documents/mfe courses/quarter2/EmpiricalMethod406/hw1/GSPC.xlsx"
DBV = as.data.frame(read_excel(file,sheet = "table (2)"))
row.names(DBV) = as.Date(DBV$Date)
DBV$Date = NULL
adj = DBV["Adj Close"]
log_return = log(adj[-1,]/adj[-length(adj[,1]),])

GSPC = as.data.frame(read_excel(gspc_file))
gspc_adj = GSPC["Adj Close"]
gspc_log_return = log(gspc_adj[-1,]/gspc_adj[-length(gspc_adj[,1]),])

dbv_mean = mean(log_return)
dbv_sigma = sd(log_return)



test_null = function(return){
  
  skewness = skewness(return)
  kurtosis = kurtosis(return)
  t_skewness = skewness / sqrt(6 * length(return))
  t_kurtosis = (kurtosis - 3)/sqrt(24/length(return))
  return(c(t_skewness,t_kurtosis))
}

dbv_skewness_t = test_null(log_return)[1]
dbv_kurtosis_t = test_null(log_return)[2]

gspc_skewness_t = test_null(gspc_log_return)[1]
gspc_kurtosis_t = test_null(gspc_log_return)[2]

test_jb = function(return){
  skewness = skewness(return)
  kurtosis = kurtosis(return)
  T = length(return)
  jb = skewness^2 / (6/T) + (kurtosis - 3)^2/(24/T)
  return(jb)
}

dbv_jb = test_jb(log_return)
gspc_jb = test_jb(gspc_log_return)


# OLS
reg_ols = lm(log_return ~ gspc_log_return)
out = summary(reg_ols)
# White standard error
white_stderr = sqrt(diag(vcovHC(reg_ols,type = "HC")))
