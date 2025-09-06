library("readxl")
market_premium = 0.06
beta = 1.4
rf = 0.05
p0 = 35
#1a
capm_return = beta * market_premium + rf

#1b
p1 = p0 * (1 + capm_return)

#1c
dividend = 2
p1_d = p0 * (1 + capm_return) - dividend

#2a
# 0.02 - rf = beta_agg * (0.05 - rf)
# 0.32 - rf = beta_agg * (0.20 - rf)

beta_agg = 2

# 0.035 - rf = beta_def * (0.05 - rf)
# 0.14 - rf = beta_def * (0.2 - rf)
beta_def = 0.7

#2b
return_agg = 0.02 * 0.5 + 0.32 * 0.5
return_def = 0.035 * 0.5 + 0.14 * 0.5

#2c
tbill = 0.08
return_agg = (tbill + beta_agg * (0.05 - tbill)) * 0.5 + 0.5 * (tbill + beta_agg * (0.2 - tbill))
return_def = (tbill + beta_def * (0.05 - tbill)) * 0.5 + 0.5 * (tbill + beta_def * (0.2 - tbill))
slope = (return_agg - return_def)/(beta_agg - beta_def)
plot(type="l",x=c(beta_agg,beta_def),y=c(return_agg,return_def),xlab = "Beta(R)",ylab = "E(R)",col = "blue",xlim = c(0,3),ylim = c(0,0.3))
abline(tbill,slope, col = "blue")

points(beta_agg,return_agg,pch=19,col="red")
points(beta_def,return_def,pch = 19, col="green")
text(beta_agg,return_agg,pch=19,col="red",labels = "Aggressive",pos = 3)
text(beta_def,return_def,pch=19,col="green",labels = "Defensive",pos = 1)
companies = c("MSFT","INTC","LUV","MCD","JNJ")

data = data.frame(msft=0,intc=0,luv=0,mcd=0,jnj=0)
data[msft,1] = 1
colnames(data) = companies
msft = read_excel("lecture8p.xlsx",sheet = "MSFT")
intc = read_excel("lecture8p.xlsx",sheet = "INTC")
luv = read_excel("lecture8p.xlsx",sheet = "LUV")
mcd = read_excel("lecture8p.xlsx", sheet = "MCD")
jnj = read_excel("lecture8p.xlsx", sheet = "JNJ")
rf = read_excel("lecture8p.xlsx", sheet = "F-F_Research_Data_Factors_daily", skip = 4)
rf = rf[rf$X__1 >= "19900101",]
rf = rf/100
data = cbind(msft$`Adj Close`,intc$`Adj Close`,luv$`Adj Close`,mcd$`Adj Close`,jnj$`Adj Close`)
colnames(data) = companies

for (i in 1:5){
  return = data[-1,i] / data[-nrow(data),i] - 1
  data[-1,i] = return
}
data = data[-1,]
msft_lm = lm((data[,1] - rf$RF) ~ rf$`Mkt-RF`)
intc_lm = lm((data[,2] - rf$RF) ~ rf$`Mkt-RF`)
luv_lm = lm((data[,3] - rf$RF) ~ rf$`Mkt-RF`)
mcd_lm = lm((data[,4] - rf$RF) ~ rf$`Mkt-RF`)
jnj_lm = lm((data[,5] - rf$RF) ~ rf$`Mkt-RF`)

# msft_beta = msft_lm$coefficients[2]
# intc_beta = intc_lm$coefficients[2]
# luv_beta = luv_lm$coefficients[2]
# mcd_beta = mcd_lm$coefficients[2]
# jnj_beta = jnj_lm$coefficients[2]


alpha_beta_msft = summary(msft_lm)$coef[,1:2]
rownames(alpha_beta_msft) = c("alpha","beta")
colnames(alpha_beta_msft) = c("estimates","standard deviation")
idiosyncratic_msft = summary(msft_lm)$sigma

alpha_beta_intc = summary(intc_lm)$coef[,1:2]
rownames(alpha_beta_intc) = c("alpha","beta")
colnames(alpha_beta_intc) = c("estimates","standard deviation")
idiosyncratic_intc = summary(intc_lm)$sigma

alpha_beta_luv = summary(luv_lm)$coef[,1:2]
rownames(alpha_beta_luv) = c("alpha","beta")
colnames(alpha_beta_luv) = c("estimates","standard deviation")
idiosyncratic_luv = summary(luv_lm)$sigma

alpha_beta_mcd = summary(mcd_lm)$coef[,1:2]
rownames(alpha_beta_mcd) = c("alpha","beta")
colnames(alpha_beta_mcd) = c("estimates","standard deviation")
idiosyncratic_mcd = summary(mcd_lm)$sigma

alpha_beta_jnj = summary(jnj_lm)$coef[,1:2]
rownames(alpha_beta_jnj) = c("alpha","beta")
colnames(alpha_beta_jnj) = c("estimates","standard deviation")
idiosyncratic_jnj = summary(jnj_lm)$sigma

highest_return = max(c(alpha_beta_msft["beta","estimates"],alpha_beta_intc["beta","estimates"],alpha_beta_luv["beta","estimates"],alpha_beta_mcd["beta","estimates"],alpha_beta_jnj["beta","estimates"]))

mean_returns = matrix(nrow = 1,ncol = 5)

mean_returns[1,] = c(alpha_beta_msft["beta","estimates"],alpha_beta_intc["beta","estimates"],alpha_beta_luv["beta","estimates"],alpha_beta_mcd["beta","estimates"],alpha_beta_jnj["beta","estimates"])

colnames(mean_returns) = companies

high_mean = names(mean_returns[,which(mean_returns[1,] == max(mean_returns[1,]))])
low_mean = names(mean_returns[,which(mean_returns[1,] == min(mean_returns[1,]))])

returns = matrix(nrow = 1,ncol = 5)
for (i in 1:5){
  returns[1,i] = mean(data[,i])
}
colnames(returns) = companies
names(returns) = companies
high = names(returns[,which(returns[1,] == max(returns[1,]))])
low = names(returns[,which(returns[1,] == min(returns[1,]))])
