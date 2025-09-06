library(DataAnalytics)
library(reshape2)
data("multi")
p1_p2 = lm(multi$p1 ~ multi$p2)

sales_e12 = lm(multi$Sales ~ p1_p2$residuals)
mr = lm(multi$Sales ~ multi$p1 + multi$p2)
sales_e12$coefficients[2] - mr$coefficients[2]


data("countryret")
y = countryret$usa
x = cbind(rep(1,length(y)),countryret$canada, countryret$uk, countryret$australia, countryret$france, countryret$germany, countryret$japan)

crossprod(x)
coefficients = chol2inv(chol(crossprod(x))) %*% crossprod(x,y)

e = y - x %*% coefficients
ssq = sum(e*e)/(length(y) - ncol(x))
Var_b = ssq * chol2inv(chol(crossprod(x)))
std_err = sqrt(diag(Var_b))
names(std_err) = c("Intercept","Canada","UK","Austrilia","France","Germany","Japan")
std_err

data("Vanguard")
data("marketRf")
Van = Vanguard[c(1,2,5)]
V_reshaped = dcast(Van,date ~ ticker, value.var = "mret")
Van_mkt = merge(V_reshaped,marketRf,by="date")
x = Van_mkt$vwretd
y = Van_mkt$VWNFX
x = x[-which(is.na(y))]
y = y[-which(is.na(y))]

out = lm(y ~ x)
s = summary(out)[[6]]
n = length(x)
s_pred = s * sqrt(1 + 1/n + (0.05 - mean(x))^2/((n - 1) * var(x)))
t = qt(0.95,df = n - 2)
interval_low = out$coefficients[1] + out$coefficients[2] * 0.05 - t * s_pred
interval_upper = out$coefficients[1] + out$coefficients[2] * 0.05 + t * s_pred


predict(out,data.frame(x = 0.05),int = "prediction",level = 0.9)



u = matrix(c(0.01,0.015,0.025),nrow = 3,ncol = 1,byrow = TRUE)
sigma = matrix(c(0.0016,0.0010,0.0015,0,0.002,0.0019,0,0,0.0042),nrow = 3,ncol = 3, byrow = TRUE)
diagnal = sqrt(diag(sigma))^-1
diagnal_matrix = diag(diagnal,nrow = 3,ncol = 3)

sigma %*% diagnal_matrix * diagnal

weights = matrix(c(0.3,0.4,0.3),nrow = 1,ncol = 3,byrow = TRUE)
weights %*% u
weights %*% sigma %*% t(weights)


mkt_up = ifelse(x>0,1,0)
mkt_down = 1 - mkt_up
mkt_up = mkt_up * x
mkt_down = mkt_down * x
mkt_timing = lm(y ~ mkt_up + mkt_down)
lmSumm(mkt_timing)
R = matrix(c(0,1,-1),byrow = TRUE,nrow = 1)
r = c(0)
x = cbind(c(rep(1,length(x))),mkt_up,mkt_down)
b = as.vector(mkt_timing$coefficients)
QFmat = chol2inv(chol(crossprod(x)))
QFmat = R %*% QFmat %*% t(R)
Violation = R%*%b - matrix(r,ncol = 1)
fnum = t(Violation) %*% chol2inv(chol(QFmat)) %*% Violation
n_minus_k = length(mkt_up) - length(b)
fdenom = nrow(R) * sum(mkt_timing$residuals ** 2)/n_minus_k
f = fnum / fdenom
f
pvalue = 1 - pf(f,df1=nrow(R),df2 = n_minus_k)
pvalue
library("quantmod")
install.packages("quantmod")
