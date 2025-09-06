library(ggplot2)
df = read.csv('company.csv')
df = na.omit(df)
gvkeys = unique(df$gvkey)
n = length(gvkeys)

make_plot = function(k){
  total_ep = rep(0, n)
  total_percent = rep(0, n)
  log_mktcap = rep(0, n)
  for (i in 1:n){
    temp = df[df$gvkey == gvkeys[i],]
    eps = temp[temp$fyear == k,]$epspi
    shares = temp[temp$fyear == k,]$csho
    price = temp[temp$fyear == k,]$prcc_f
    mktcap = log(temp[temp$fyear == k,]$csho * temp[temp$fyear == k,]$prcc_f)
    EP_ratio = eps/price
    earning_k = temp[temp$fyear == k,]$ebitda
    earning_k_3 = temp[temp$fyear == (k-3), ]$ebitda
    percent = earning_k/earning_k_3 - 1
    if (length(EP_ratio)>0 & length(percent)>0 & length(mktcap) > 0){
      total_ep[i] = EP_ratio
      total_percent[i] = percent
      log_mktcap[i] = mktcap
    }else{
      total_ep[i] = NA
      total_percent[i] = NA
      log_mktcap[i] = NA
    }
  }
  out = as.data.frame(cbind(total_ep,total_percent, log_mktcap))
  out = na.omit(out)
  int = quantile(out$total_percent, c(0.1,0.9))
  int2 = quantile(out$total_ep, c(0.1,0.9))
  outsub = out[out$total_percent > int[[1]] & out$total_percent < int[[2]],]
  outsub = outsub[outsub$total_ep > int2[[1]] & outsub$total_ep < int2[[2]],]
  outsub = outsub[outsub$log_mktcap > 0,]
  
  ggplot(outsub, aes(total_percent, total_ep)) + geom_point(aes(size = log_mktcap)) + 
    ggtitle(paste('E/P ratio to earnings change percentage of each company for fiscal year ',k)) +
    xlab(paste('earnings change percentage between fiscal year ',k,' and ', k-3)) +
    ylab('E/P ratio')
}

make_plot(2018)
make_plot(2017)
make_plot(2016)
make_plot(2015)
