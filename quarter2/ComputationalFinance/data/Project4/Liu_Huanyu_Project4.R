df = read.csv("GOOGL.csv")
price = df[[6]]
return = price[-1] / price[-length(price)] - 1
std = sd(return)
std = sqrt(252) * std