eth = read.csv("eth2.csv")
eth$date = as.Date(as.character(eth$date), format = "%Y%m%d")
market_rf = read.csv("F-F_Research_Data_Factors_daily.CSV",skip = 4)
market_rf$X = as.Date(as.character(market_rf$X),format = "%Y%m%d")
market_rf = market_rf[-nrow(market_rf),]
rownames(eth) = eth$date
rownames(market_rf) = market_rf$X
eth$date = NULL
market_rf$X = NULL
eth = eth[rownames(eth) >= "2016-09-30" & rownames(eth) <= "2017-09-29",]
market_rf = market_rf/100
mean_mkt_rf = mean(market_rf$Mkt.RF)
mean_smb = mean(market_rf$SMB)
mean_hml = mean(market_rf$HML)
market_rf = market_rf[rownames(market_rf) >= "2016-09-30" & rownames(market_rf) <= "2017-09-29",]

risk_free = 2.33 / 252 /100
return_premium = as.numeric(levels(eth$RET))[eth$RET] - market_rf$RF
capm = lm(return_premium ~ market_rf$Mkt.RF)
risk_free = 2.33 / 252 /100
capm_return = capm$coefficients[2] * mean_mkt_rf + risk_free
annual_capm = (1 + capm_return) ^ 252 - 1
factor_3 = lm(return_premium ~ market_rf$Mkt.RF + market_rf$SMB + market_rf$HML)
factor_3_return = factor_3$coefficients[2] * mean_mkt_rf + factor_3$coefficients[3] * mean_smb + factor_3$coefficients[4] * mean_hml + risk_free
annual_3_factor = (1 + factor_3_return) ^ 252 - 1
