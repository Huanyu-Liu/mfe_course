import pandas as pd
import numpy as np
SMALL = 0.00001

zcb = pd.read_excel('/Users/huanyu/Desktop/FixedIncome/hw1/zcb.xlsx')
cb = pd.read_excel('/Users/huanyu/Desktop/FixedIncome/hw1/cb.xlsx')
cb = pd.merge(cb,zcb,how='left',on='Maturity')
cb.Maturity.duplicated()
selected_maturity = np.sort(np.concatenate((cb[cb.Maturity.duplicated()].index, cb[cb.Maturity.duplicated()].index - 1)))
cb = cb.loc[selected_maturity,:]
result = list()
for i in range(0,len(cb),2):
    x = 1
    y = -cb.iloc[i,1] * x / cb.iloc[i+1,1]
    z = -y - x
    synthetic_price = x * cb.iloc[i,2] + y * cb.iloc[i+1,2]
    zcb_price = z * cb.iloc[i,3]
    initial_cost = synthetic_price + zcb_price
    if (abs(initial_cost) > SMALL):
        if initial_cost > 0:
            result.append([-x,cb.iloc[i,1], -y,cb.iloc[i+1,1], -z, cb.iloc[i, 0],initial_cost])
        else:
            result.append([x,cb.iloc[i,1], y,cb.iloc[i+1,1], z, cb.iloc[i, 0],initial_cost])
result_df = pd.DataFrame({'Coupon Bond A':[x[0] for x in result],'Coupon A':[x[1] for x in result],'Coupon Bond B':[x[2] for x in result],'Coupon B':[x[3] for x in result],'Zero Coupon Bond':[x[4] for x in result],'Maturity':[x[5] for x in result],'Arbitrage Profit':[x[6] for x in result]})