# Problem 1

industries = read.csv('48_Industry_Portfolios.CSV')
ff_research = read.csv('F-F_Research_Data_Factors.CSV')
industries[industries == -99.99] = NA
industries = industries[ , colSums(is.na(industries)) == 0]
excess_return = sweep(industries[,-1], 1, ff_research[,'RF'], '-')
varcov = cov(excess_return)
eigen_value = eigen(varcov)
sum_eigen = sum(eigen_value$values)
frac = eigen_value$values / sum_eigen
barplot(frac)

# Problem 2
# (a)
sorted = sort(frac, decreasing = T)
total3_largest = sorted[1] + sorted[2] + sorted[3]
total3_largest

# (b)
#summary(princomp(excess_return))
prcomp_result = prcomp(excess_return)
weight_3_factor = prcomp_result$rotation[,1:3]
sd_3_factor = prcomp_result$sdev[1:3]
returns_3_factor = as.matrix(excess_return) %*% as.matrix(weight_3_factor)
mean_returns = apply(returns_3_factor, 2, FUN = mean)
cor_3_factor = cor(returns_3_factor)


# (c)
betas = eigen_value$vectors[,c(1:3)]
expected_return = as.matrix(excess_return) %*% as.matrix(betas)
predict_return = as.matrix(expected_return) %*% as.matrix(t(betas))
predict_return_mean = apply(predict_return, 2, mean)
sample_industry_mean = apply(excess_return, 2, mean)
plot(x = sample_industry_mean, y = predict_return_mean, ylim = c(0,1.2), xlim = c(0,1.2), pch = 19, col = 'red', cex = 0.6)
abline(0,1)

# (d)
R_sqr = 1 - var(sample_industry_mean - predict_return_mean) / var(sample_industry_mean)
R_sqr
# Problem 3
# (a)
industries25 = read.csv('25_Portfolios_5x5.csv', header = T)
excess_return25 = industries25[2:26] - ff_research$RF
varcov25 = cov(excess_return25)
eigen_value25 = eigen(varcov25)
sum_eigen25 = sum(eigen_value25$values)
frac25 = eigen_value25$values / sum_eigen25
barplot(frac25)

# (b)

summary(princomp(excess_return25))

# I would say 4 factors should explain average returns to the 25 F-F portfolios.