require(data.table)
require(ggplot2)
require(foreign)
lending_club = as.data.table(read.dta('/Users/huanyu/Desktop/MachineLearning431/hw3/LendingClub_LoanStats3a_v12.dta'))
lending_club = lending_club[lending_club$loan_status == 'Fully Paid' | lending_club$loan_status == 'Charged Off']
lending_club[lending_club$loan_status == 'Fully Paid','Default'] = 0
lending_club[lending_club$loan_status == 'Charged Off','Default'] = 1
average_default_rate = length(which(lending_club$Default == 1)) / length(lending_club$Default)
regression_output = glm(Default~grade,data=lending_club,family="binomial")
summary(regression_output)
