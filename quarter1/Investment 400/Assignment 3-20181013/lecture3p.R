#Problem Set 3 Solution Template
# mysoln is a list to store your answers

# the first element of the list is a vector of names for the students in your group
# make sure these match the names shown on the pdf document provided by the MFE office
# using group 1 as an example:
mysoln = list(student = c("Molin Liang", "Meghana Rao", "Chengbo Du", "Shardul Kulkarni"))



# 4

# your intermediary code to get your answers here

bondA.weight = 2000000 / 100

A = matrix(data=c(4, 6, 104, 106), nrow=2, ncol=2, byrow=TRUE)    
b = matrix(data=c(0,1000000),nrow=2, ncol=1, byrow=FALSE)
b.and.c = solve(A, b)
bondB.weight = b.and.c[1]
bondC.weight = b.and.c[2]
# answers
# Note: Remember to add conclusions in write-up

discount = 1 + 0.06
price.D = 100 / discount^10
price.E = 100 / discount^15
liability = 1000000 / discount^30 + 2000000 / discount^31
weight.D = price.D / liability
weight.E = price.E / liability
liability_duration = 1000000/discount^30 / liability * 30 + 2000000/discount^31 / liability * 31
A = matrix(data = c(10 * weight.D, 15 * weight.E, price.D, price.E), nrow = 2, ncol = 2, byrow = TRUE)
b = matrix(data = c(liability_duration, liability), nrow = 2, ncol = 1, byrow = FALSE)
D.and.E = solve(A,b)
bondD.weight = D.and.E[1]
bondE.weight = D.and.E[2]



a = c(bondA.weight, bondB.weight, bondC.weight)
b = c(bondD.weight, bondE.weight)


#c = "Put in PDF writeup"
#d = "Put in PDF writeup"

# add answers to list for "Q4"
mysoln[["Q4"]] =  list(a=a, b=b) #c and d in writeup


# return my solution
mysoln