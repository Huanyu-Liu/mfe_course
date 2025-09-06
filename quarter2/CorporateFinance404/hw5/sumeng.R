library(quantmod)
library(lubridate)
library(dplyr)
library(ggplot2)

df = read.csv("hw3.csv")
perm_num = unique(df$PERMNO)
n = length(perm_num)
result = list()
ret_2008 = c()
ret_2009 = c()
ret_2010 = c()
ret_2011 = c()
ret_2012 = c()
ret_2013 = c()
ret_2014 = c()
ret_2015 = c()
ret_2016 = c()
ret_2017 = c()
ret_2018 = c()
mktcap_2008 = c()
mktcap_2009 = c()
mktcap_2010 = c()
mktcap_2011 = c()
mktcap_2012 = c()
mktcap_2013 = c()
mktcap_2014 = c()
mktcap_2015 = c()
mktcap_2016 = c()
mktcap_2017 = c()
mktcap_2018 = c()
for (i in 1:n){
  df_1 = df[df$PERMNO == perm_num[i],]
  df_xts = as.xts(df_1, order.by = as.Date(as.character(df_1$date), "%Y%m%d"))
  t = df_xts[endpoints(df_xts, 'years'),]
  ret = diff(as.numeric(t$PRC)) / as.numeric(t[1:length(t$PERMNO)-1]$PRC)
  yr = year(t[-1]$date)
  out = as.data.frame(t[-1])
  mktcap = as.numeric(t[-1]$PRC) * as.numeric(t[-1]$SHROUT)
  out$RET = ret
  out$Year = yr
  out$mktcap = mktcap
  
  ret_2008 = if (length(out[out$Year == 2008, "RET"])>0) c(ret_2008, out[out$Year == 2008, "RET"]) else c(ret_2008, NA)
  ret_2009 = if (length(out[out$Year == 2009, "RET"])>0) c(ret_2009, out[out$Year == 2009, "RET"]) else c(ret_2009, NA)
  ret_2010 = if (length(out[out$Year == 2010, "RET"])>0) c(ret_2010, out[out$Year == 2010, "RET"]) else c(ret_2010, NA)
  ret_2011 = if (length(out[out$Year == 2011, "RET"])>0) c(ret_2011, out[out$Year == 2011, "RET"]) else c(ret_2011, NA)
  ret_2012 = if (length(out[out$Year == 2012, "RET"])>0) c(ret_2012, out[out$Year == 2012, "RET"]) else c(ret_2012, NA)
  ret_2013 = if (length(out[out$Year == 2013, "RET"])>0) c(ret_2013, out[out$Year == 2013, "RET"]) else c(ret_2013, NA)
  ret_2014 = if (length(out[out$Year == 2014, "RET"])>0) c(ret_2014, out[out$Year == 2014, "RET"]) else c(ret_2014, NA)
  ret_2015 = if (length(out[out$Year == 2015, "RET"])>0) c(ret_2015, out[out$Year == 2015, "RET"]) else c(ret_2015, NA)
  ret_2016 = if (length(out[out$Year == 2016, "RET"])>0) c(ret_2016, out[out$Year == 2016, "RET"]) else c(ret_2016, NA)
  ret_2017 = if (length(out[out$Year == 2017, "RET"])>0) c(ret_2017, out[out$Year == 2017, "RET"]) else c(ret_2017, NA)
  ret_2018 = if (length(out[out$Year == 2018, "RET"])>0) c(ret_2018, out[out$Year == 2018, "RET"]) else c(ret_2018, NA)
  
  mktcap_2008 = if (length(out[out$Year == 2008, "mktcap"])>0) c(mktcap_2008, out[out$Year == 2008, "mktcap"]) else c(mktcap_2008, NA)
  mktcap_2009 = if (length(out[out$Year == 2009, "mktcap"])>0) c(mktcap_2009, out[out$Year == 2009, "mktcap"]) else c(mktcap_2009, NA)
  mktcap_2010 = if (length(out[out$Year == 2010, "mktcap"])>0) c(mktcap_2010, out[out$Year == 2010, "mktcap"]) else c(mktcap_2010, NA)
  mktcap_2011 = if (length(out[out$Year == 2011, "mktcap"])>0) c(mktcap_2011, out[out$Year == 2011, "mktcap"]) else c(mktcap_2011, NA)
  mktcap_2012 = if (length(out[out$Year == 2012, "mktcap"])>0) c(mktcap_2012, out[out$Year == 2012, "mktcap"]) else c(mktcap_2012, NA)
  mktcap_2013 = if (length(out[out$Year == 2013, "mktcap"])>0) c(mktcap_2013, out[out$Year == 2013, "mktcap"]) else c(mktcap_2013, NA)
  mktcap_2014 = if (length(out[out$Year == 2014, "mktcap"])>0) c(mktcap_2014, out[out$Year == 2014, "mktcap"]) else c(mktcap_2014, NA)
  mktcap_2015 = if (length(out[out$Year == 2015, "mktcap"])>0) c(mktcap_2015, out[out$Year == 2015, "mktcap"]) else c(mktcap_2015, NA)
  mktcap_2016 = if (length(out[out$Year == 2016, "mktcap"])>0) c(mktcap_2016, out[out$Year == 2016, "mktcap"]) else c(mktcap_2016, NA)
  mktcap_2017 = if (length(out[out$Year == 2017, "mktcap"])>0) c(mktcap_2017, out[out$Year == 2017, "mktcap"]) else c(mktcap_2017, NA)
  mktcap_2018 = if (length(out[out$Year == 2018, "mktcap"])>0) c(mktcap_2018, out[out$Year == 2018, "mktcap"]) else c(mktcap_2018, NA)
  result[[i]] = out
}

mean_ret_2008 = mean(ret_2008, na.rm = T)
mean_ret_2009 = mean(ret_2009, na.rm = T)
mean_ret_2010 = mean(ret_2010, na.rm = T)
mean_ret_2011 = mean(ret_2011, na.rm = T)
mean_ret_2012 = mean(ret_2012, na.rm = T)
mean_ret_2013 = mean(ret_2013, na.rm = T)
mean_ret_2014 = mean(ret_2014, na.rm = T)
mean_ret_2015 = mean(ret_2015, na.rm = T)
mean_ret_2016 = mean(ret_2016, na.rm = T)
mean_ret_2017 = mean(ret_2017, na.rm = T)
mean_ret_2018 = mean(ret_2018, na.rm = T)

mean_mktcap_2008 = mean(mktcap_2008, na.rm = T)
mean_mktcap_2009 = mean(mktcap_2009, na.rm = T)
mean_mktcap_2010 = mean(mktcap_2010, na.rm = T)
mean_mktcap_2011 = mean(mktcap_2011, na.rm = T)
mean_mktcap_2012 = mean(mktcap_2012, na.rm = T)
mean_mktcap_2013 = mean(mktcap_2013, na.rm = T)
mean_mktcap_2014 = mean(mktcap_2014, na.rm = T)
mean_mktcap_2015 = mean(mktcap_2015, na.rm = T)
mean_mktcap_2016 = mean(mktcap_2016, na.rm = T)
mean_mktcap_2017 = mean(mktcap_2017, na.rm = T)
mean_mktcap_2018 = mean(mktcap_2018, na.rm = T)

sd_ret_2008 = sd(ret_2008, na.rm = T)
sd_ret_2009 = sd(ret_2009, na.rm = T)
sd_ret_2010 = sd(ret_2010, na.rm = T)
sd_ret_2011 = sd(ret_2011, na.rm = T)
sd_ret_2012 = sd(ret_2012, na.rm = T)
sd_ret_2013 = sd(ret_2013, na.rm = T)
sd_ret_2014 = sd(ret_2014, na.rm = T)
sd_ret_2015 = sd(ret_2015, na.rm = T)
sd_ret_2016 = sd(ret_2016, na.rm = T)
sd_ret_2017 = sd(ret_2017, na.rm = T)
sd_ret_2018 = sd(ret_2018, na.rm = T)

sd_mktcap_2008 = sd(mktcap_2008, na.rm = T)
sd_mktcap_2009 = sd(mktcap_2009, na.rm = T)
sd_mktcap_2010 = sd(mktcap_2010, na.rm = T)
sd_mktcap_2011 = sd(mktcap_2011, na.rm = T)
sd_mktcap_2012 = sd(mktcap_2012, na.rm = T)
sd_mktcap_2013 = sd(mktcap_2013, na.rm = T)
sd_mktcap_2014 = sd(mktcap_2014, na.rm = T)
sd_mktcap_2015 = sd(mktcap_2015, na.rm = T)
sd_mktcap_2016 = sd(mktcap_2016, na.rm = T)
sd_mktcap_2017 = sd(mktcap_2017, na.rm = T)
sd_mktcap_2018 = sd(mktcap_2018, na.rm = T)

ret_to_cap = cbind(c(mean_ret_2008,
                     mean_ret_2009,
                     mean_ret_2010,
                     mean_ret_2011,
                     mean_ret_2012,
                     mean_ret_2013,
                     mean_ret_2014,
                     mean_ret_2015,
                     mean_ret_2016,
                     mean_ret_2017,
                     mean_ret_2018),
                   c(mean_mktcap_2008,
                     mean_mktcap_2009,
                     mean_mktcap_2010,
                     mean_mktcap_2011,
                     mean_mktcap_2012,
                     mean_mktcap_2013,
                     mean_mktcap_2014,
                     mean_mktcap_2015,
                     mean_mktcap_2016,
                     mean_mktcap_2017,
                     mean_mktcap_2018))
ret_to_cap = as.data.frame(ret_to_cap)
names(ret_to_cap) = c('ret','mktcap')
ggplot(data = ret_to_cap, aes(y = ret, x = mktcap)) + geom_point() + ggtitle('average return to average market cap of all stocks')