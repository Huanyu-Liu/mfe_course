library(foreign)
library(data.table)
df = read.dta('LendingClub_LoanStats3a_v12.dta')
# a
# (i)
df = df[df$loan_status == 'Fully Paid' | df$loan_status == 'Charged Off',]
df$default = 0
df[df$loan_status == 'Charged Off', 'default'] = 1
# (ii)
num_defaults = sum(df$default)
default_rate = num_defaults / nrow(df)

# b
# (i)
out = glm(formula = default ~ grade, family = 'binomial', data = df)
summary(out)
# (ii)
test_stat = out$null.deviance - out$deviance
k = out$df.null - out$df.residual
pval_chisq = 1 - pchisq(test_stat,df = k)
# (iii)
phat = predict(out, type = 'response')
ranks = rank(phat, ties.method = 'first')
deciles = cut(ranks, quantile(ranks, probs=0:10/10), include.lowest=T)
deciles = as.numeric(deciles)
df2 = data.frame(deciles=deciles,phat=phat,default=df$default)
lift=aggregate(df2,by=list(deciles),FUN="mean",data=df2)
lift=lift[,c(2,4)]
lift[,3]=lift[,2]/mean(df$default)
names(lift)=c("decile","Mean Response","Lift Factor")
lift

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}
roc = simple_roc(df$default == '1', phat)
TPR = roc$TPR
FPR = roc$FPR
library(ggplot2)
q = qplot(FPR,TPR,xlab = 'FPR', ylab = 'TPR', col = I('red'), main = 'ROC Curve, Logistic Default Model', size = I(0.5))
q = q + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1)) + theme_bw()
q
#(iv)
profits = sum(df$default == '0')*(1-FPR)*1
loss = sum(df$default == '1')*(1-TPR)*10
net = profits - loss
index = which.max(net)
cutoff = phat[index]
q + geom_point(aes(x = FPR[index], y = TPR[index], col = I('green'), size = I(2)))
#c
#(i)
out2 = glm(formula = default ~ loan_amnt + annual_inc, family = 'binomial', data = df)
summary(out2)
phat2 = predict(out2, type = 'response')
deciles2 = cut(phat2, quantile(phat2, probs=0:10/10), include.lowest=T)
deciles2 = as.numeric(deciles2)
df3 = data.frame(deciles=deciles2,phat=phat2,default=df$default)
lift2=aggregate(df3,by=list(deciles2),FUN="mean",data=df2)
lift2=lift2[,c(2,4)]
lift2[,3]=lift2[,2]/mean(df$default)
names(lift2)=c("decile","Mean Response","Lift Factor")
lift2
roc2 = simple_roc(df$default == '1', phat2)
TPR2 = roc2$TPR
FPR2 = roc2$FPR
roc <- cbind(roc,Model = "V1")
roc2 <- cbind(roc2, Model = "V2")

New_ROC <- rbind(roc, roc2)
q = qplot(FPR,TPR, data = New_ROC, xlab = 'FPR', ylab = 'TPR', col = Model, main = 'ROC Curve, Logistic Default Model 2', size = I(0.5))
q + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), col = I('black')) + theme_bw()
#(ii)
out3 = glm(formula = default ~ loan_amnt + annual_inc + term + int_rate, family = 'binomial', data = df)
summary(out3)
phat3 = predict(out3, type = 'response')
deciles3 = cut(phat3, quantile(phat3, probs=0:10/10), include.lowest=T)
deciles3 = as.numeric(deciles3)
df4 = data.frame(deciles=deciles3,phat=phat3,default=df$default)
lift3=aggregate(df4,by=list(deciles3),FUN="mean",data=df4)
lift3=lift3[,c(2,4)]
lift3[,3]=lift3[,2]/mean(df$default)
names(lift3)=c("decile","Mean Response","Lift Factor")
lift3
roc3 = simple_roc(df$default == '1', phat3)
TPR3 = roc3$TPR
FPR3 = roc3$FPR
roc3 <- cbind(roc3, Model = "V3")

New_ROC <- rbind(New_ROC, roc3)
q = qplot(FPR,TPR, data = New_ROC, xlab = 'FPR', ylab = 'TPR', col = Model, main = 'ROC Curve, Logistic Default Model 2', size = I(0.5))
q + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), col = I('black')) + theme_bw()

#(iii)
df$int_rate_sq = df$int_rate ^ 2
out4 = glm(formula = default ~ loan_amnt + annual_inc + term + int_rate + int_rate_sq, family = 'binomial', data = df)
summary(out4)

