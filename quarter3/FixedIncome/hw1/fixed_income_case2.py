import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
data = pd.read_excel('/Users/huanyu/Desktop/FixedIncome/hw1/Data for Case 2.xlsx', skiprows=5, header=None, usecols=[1,2,3])
data.columns = ['date', 'principal_strip', 'coupon_strip']
trading_days = len(data)
principal_unit = 100 / data.loc[0,'principal_strip']
coupon_unit = -100 / data.loc[0,'coupon_strip']

# long principal strip and short coupon strip
r = 0.002
data['portfolio_value'] = principal_unit * data['principal_strip'] + coupon_unit * data['coupon_strip'] + np.array([5 * (1+r/trading_days)**x for x in range(0,252)])

data['return'] = np.divide(data.loc[1:,'portfolio_value'],data.iloc[:-1,3]) - 1
returns = data.loc[1:,'return']
plt.plot(data['portfolio_value'])
plt.show()
distribution = sns.distplot(returns, bins=20)
plt.show()
realized_return = data.iloc[-1,3] / 5 - 1
volatility = returns.std()
annual_v = volatility * math.sqrt(trading_days)
Annual_SR = (realized_return - r) / annual_v
excess_kurtosis = returns.kurtosis()
max_negative_ret = returns.min()
