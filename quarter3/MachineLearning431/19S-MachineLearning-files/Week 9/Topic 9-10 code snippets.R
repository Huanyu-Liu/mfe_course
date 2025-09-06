#
# code snippets for Topic 9-10
#
# set your working directory, where you have downloaded the Stata data file
# change the below to match your folder
#setwd("D:/llochsto/Dropbox/Data Analytics/Data")
setwd("C:/Users/Lars/Dropbox (Personal)/Data Analytics/Data")

require(data.table)
require(ggplot2)


#
# let's look at the loan data again
#
require(DataAnalytics)
data(loans)
loans$default=loans$Good.Loan-1
loans$Good.Loan=NULL # remove old variable
names(loans)=c("StatChkA","Duration","CrdHist","Purpose","CrdAmt","Sav_Bnd",
               "Emply","InstallRate","Pstatus","OthrDebt",
               "Resid","Proprty","Age","OthrInstall","Housing","Ncredits","job","Nsupport",
               "Telephone","Foreign","default")

# let's stripdown the names so that we can view the output on one "page"

#
# three types of variables
# Credit History
# loan characteristics
# demographics
#
# let's fit logistic regression with only the factors that have significant coefficients 
#
x_data = loans[, c("default","StatChkA","CrdHist","Sav_Bnd","OthrDebt","OthrInstall","Purpose","Duration","CrdAmt","InstallRate","Pstatus","Foreign")]
outloans1=glm(default~.,data=x_data,family="binomial")

#
# lets look at the fitted probabilities 
#
phat=predict(outloans1,type="response")

# compute a ROC curve
# define a function that creates the true and false positive rates
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

glm_simple_roc <- simple_roc(loans$default=="1", phat)
TPR <- glm_simple_roc$TPR
FPR <- glm_simple_roc$FPR

# plot the corresponding ROC curve
q <- qplot(FPR,TPR,xlab="FPR",ylab="TPR",col=I("blue"),main="ROC Curve for Logistic Regression Default Model",size=I(0.75))
# add straight 45 degree line from 0 to 1
q + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size=I(1.0)) + theme_bw()



# now, fit an SVM to the same data and compare ROC curves
require(e1071)

# first, do a polynomial kernel
outloans_svm_poly = svm(default~.,data=x_data,kernel = "polynomial", degree = 3, gamma = 1/11, probability = TRUE)
phat_svm_poly=predict(outloans_svm_poly,x_data, probability = TRUE)

svm_poly_roc <- simple_roc(x_data$default=="1", phat_svm_poly)
TPR <- svm_poly_roc$TPR
FPR <- svm_poly_roc$FPR

# plot the corresponding ROC curve
q_poly <- qplot(FPR,TPR,xlab="FPR",ylab="TPR",col=I("blue"),main="ROC Curve for Polynomial SVM",size=I(0.75))
# add straight 45 degree line from 0 to 1
q_poly + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size=I(1.0)) + theme_bw()

# second, do radial SVM
outloans_svm_radi = svm(default~.,data=x_data,kernel = "radial", gamma = 1/11, probability = TRUE)
phat_svm_radi=predict(outloans_svm_radi,x_data, probability = TRUE)

svm_radi_roc <- simple_roc(x_data$default=="1", phat_svm_radi)
TPR <- svm_radi_roc$TPR
FPR <- svm_radi_roc$FPR

# plot the corresponding ROC curve
q_radi <- qplot(FPR,TPR,xlab="FPR",ylab="TPR",col=I("blue"),main="ROC Curve for Radial SVM",size=I(0.75))
# add straight 45 degree line from 0 to 1
q_radi + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size=I(1.0)) + theme_bw()


