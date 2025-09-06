#1
#Sample Mean
Avg.Mkt = 0.05
Avg.1 = 0.01 + 0.9*Avg.Mkt
Avg.2 = -0.015 + 1.2*Avg.Mkt
Avg.3 = 0.005 + 1.0*Avg.Mkt
Avg = c(Avg.1,Avg.2,Avg.3)
Avg_df = data.frame(Avg.1,Avg.2,Avg.3)
colnames(Avg_df) = c("1st","2nd","3rd")
Avg_df
#Standard Deviation
Variance_matrix = matrix(c(0.1^2,0,0,0,0.15^2,0,0,0,0.05^2),nrow=3)
betas = c(0.9,1.2,1)
Variance.Mkt = 0.15^2
Variance.1 = betas[1]^2 * Variance.Mkt +  Variance_matrix[1,1] 
Variance.2 = betas[2]^2 * Variance.Mkt +  Variance_matrix[2,2]
Variance.3 = betas[3]^2 * Variance.Mkt +  Variance_matrix[3,3]
Sd = c(sqrt(Variance.1),sqrt(Variance.2),sqrt(Variance.3))
Sd_df = data.frame(sqrt(Variance.1),sqrt(Variance.2),sqrt(Variance.3))
colnames(Sd_df) = c("1st","2nd","3rd")
Sd_df
#SharpeRatio
Avg/Sd
#2
#Sample Mean
Avg.Hedged.1 = 0.01
Avg.Hedged.2 = -0.015
Avg.Hedged.3 = 0.005
Avg.Hedged = c(Avg.Hedged.1,Avg.Hedged.2,Avg.Hedged.3)
Avg.Hedged_df = data.frame(Avg.Hedged.1,Avg.Hedged.2,Avg.Hedged.3)
colnames(Avg.Hedged_df) = c("1st Hedged","2nd Hedged","3rd Hedged")
Avg.Hedged_df
#Standard Deviation
Variance.Hedged.1 = Variance_matrix[1,1] 
Variance.Hedged.2 = Variance_matrix[2,2]
Variance.Hedged.3 = Variance_matrix[3,3]
Sd.Hedged = c(sqrt(Variance.Hedged.1),sqrt(Variance.Hedged.2),sqrt(Variance.Hedged.3))
Sd.Hedged_df = data.frame(sqrt(Variance.Hedged.1),sqrt(Variance.Hedged.2),sqrt(Variance.Hedged.3))
colnames(Sd.Hedged_df) = c("1st Hedged","2nd Hedged","3rd Hedged")
Sd.Hedged_df
#SharpeRatio
Avg.Hedged/Sd.Hedged
#3
SharpeRatioSq.Max = t(Avg.Hedged)%*%chol2inv(chol(Variance_matrix))%*%Avg.Hedged
sqrt(SharpeRatioSq.Max)
#4 
SharpeRatio.Market = 1/3
SharpeRatio.Combined = sqrt(SharpeRatioSq.Max + SharpeRatio.Market^2)
SharpeRatio.Combined
#5
##5a
AllReturns = c(Avg,Avg.Mkt)
#Calculate systematic variance and covariance (Betai * Betaj * market variance)
betas5 = c(betas,1)
systematicVar.5 = (betas5%*%t(betas5))*Variance.Mkt
fullCovarianceMatrix.5 = rbind(cbind(Variance_matrix,0),0) + systematicVar.5
SharpeRatio.Max.Combined = t(AllReturns)%*%chol2inv(chol(fullCovarianceMatrix.5))%*%AllReturns
Sd = 0.15
k = Sd/sqrt(SharpeRatio.Max.Combined)
weights_combined = (chol2inv(chol(fullCovarianceMatrix.5))%*%AllReturns) * as.numeric(k)
rownames(weights_combined) = c("Stock1","Stock2","Stock3","Market")
weights_combined
##5b  
#Mean
Mean5 = AllReturns%*%weights_combined
#SD
Sd5 = sqrt(t(weights_combined)%*%fullCovarianceMatrix.5%*%weights_combined)
#Sharpe Ratio
SR5 = Mean5/Sd5
output = data.frame(Mean5,Sd5,SR5)
names(output) = c("Mean","SD","Sharpe Ratio")
output
#6
##6a
mimick.Weights = ((betas - mean(betas))/(length(betas)*(mean(betas^2)-mean(betas)^2)))
mimick.Return = mimick.Weights%*%Avg
systematicVar.stocks = (betas%*%t(betas))*Variance.Mkt
fullCovarianceMatrix.stocks = systematicVar.stocks +Variance_matrix
mimick.Sd = sqrt(t(mimick.Weights)%*%fullCovarianceMatrix.stocks%*%mimick.Weights)
mimick.sharpe = mimick.Return/mimick.Sd
output = data.frame(mimick.Return,mimick.Sd,mimick.sharpe)
names(output) = c("Mean","SD","Sharpe Ratio")
output
##6b
cor.mimick.market = mimick.Weights%*%(betas*sqrt(Variance.Mkt)/as.vector(mimick.Sd))
cor.mimick.market
##6c
eigens = eigen(fullCovarianceMatrix.stocks)
output = eigens$values/sum(eigens$values)
names(output) = c("1st PCA","2nd PCA","3rd PCA")
output
##6d
#Weights
portfolioweights = eigens$vectors
colnames(portfolioweights) = c("1st PCA","2nd PCA","3rd PCA")
row.names(portfolioweights) = c("1st Stock","2nd Stock","3rd Stock")
portfolioweights
#Loadings
loadings = matrix(nrow=3,ncol=3)
for(i in 1:length(eigens$values)){
    loadings[,i] = portfolioweights[,i]*sqrt(eigens$values[i])
}
colnames(loadings) = c("1st PCA","2nd PCA","3rd PCA")
row.names(loadings) = c("1st Stock","2nd Stock","3rd Stock")
loadings