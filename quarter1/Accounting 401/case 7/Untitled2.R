eth = read.csv("eth2.csv")
market_rf = read.csv("F-F_Research_Data_Factors_daily.CSV",skip = 4)
eth$date = as.Date(as.character(eth$date), "%Y%m%d")
rownames(eth) = eth$date
eth$date = NULL
eth = eth[-1,]
market_rf = market_rf[-length(market_rf[,1]),]
market_rf$X = as.Date(as.character(market_rf$X),"%Y%m%d")

rownames(market_rf) = market_rf$X
market_rf = market_rf[market_rf$X >= "1993-03-17" & market_rf$X <= "2018-09-28",]
market_rf$Mkt.RF = market_rf$Mkt.RF /100
a = market_rf$Mkt.RF
b = as.numeric(levels(eth$RET))[eth$RET]
result = lm(b ~ a)
summary(result)
result2 = lm(b ~ market_rf$Mkt.RF + market_rf$SMB + market_rf$HML)
