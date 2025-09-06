# Problem Set 2 Solution Template
# mysoln is a list to store your answers

# the first element of the list is a vector of names for the students in your group
# make sure these match the names shown on the pdf document provided by the MFE office
# using group 1 as an example:
mysoln = list(student = c("Molin Liang", "Meghana Rao", "Chengbo Du", "Shardul Kulkarni"))
# 1

# your intermediary code to get your answers here

# save down your final answers for part a, b, and c
price_func = function(coupon_rate, year, ytm){
  price = 0
  for (i in 1:year){
    price = price + 100 * coupon_rate / (1 + ytm) ^ i
  }
  price = price + 100 / (1 + ytm) ^ year
  return(price)
}
ytm = 0.06
coupon_rate = 0
year = 3
price = price_func(coupon_rate = 0, year = 3, ytm = 0.06)
a = c(price,ytm) #ytm in decimal form

ytm = 0.055
coupon_rate = 0.06
year = 2
price = price_func(coupon_rate = 0.06, year = 2, ytm = 0.055)
b = c(price,ytm) #ytm in decimal form

ytm = 0.063
coupon_rate = 0.08
year = 4
price = price_func(coupon_rate = 0.08, year = 4, ytm = 0.063)
c = c(price,ytm) #ytm in decimal form

# add answers to list for "Q1"
mysoln[["Q1"]] = list(a=a, b=b, c=c)

# 2

# your intermediary code to get your answers here

# save down your final answers

couponx = 100 * 0.04 * 0.5
pricex = 100.98
par = 100
r.6month = ((par + couponx) / pricex - 1) * 2
coupony = 100 * 0.06 * 0.5
pricey = 103.59

r.1yr = (par + coupony) / (pricey - coupony / (1 + r.6month)) - 1
a = c(r.6month, r.1yr) #in decimal form

# add answers to list for "Q2"
mysoln[["Q2"]] = list(a=a)

# 3
# your intermediary code to get your answers here

# save down your final answers
# suppose we buy x Bond A, y Bond B, z Bond C
# To eliminate any future payments, we have
# 100x + 5y + 7z = 0
# 105y + 107z = 0
# x = (-2/105)z
# y = (-105/107)z
# suppose we short 1 Bond C
priceA = 95.238
priceB = 98.438
priceC = 103.370

arbitrage = -2/105 * priceA - 105/107 * priceB + priceC
arbitrage
# Put the answer in your PDF writeup
mysoln
