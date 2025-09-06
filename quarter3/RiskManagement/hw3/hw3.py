import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.stats import norm
from scipy.stats import chi2
import seaborn as sns
# 1
c = 0.99
data = pd.read_csv('/Users/huanyu/Desktop/RiskManagement/hw3/hw3_returns2.csv')
data['Date'] = pd.to_datetime(data['Date'])
data.loc[252:,'history_VaR'] = data.index[252:].map(lambda x: np.sort(data.loc[:x-1,'Return'])[math.ceil(x * 0.01) - 1])
lamb = 0.995
length = len(data)
for i in range(252,length):
    weight = [lamb ** (i - j) * (1 - lamb) / (1 - lamb ** i) for j in range(1,i+1)]
    sorted_data = pd.DataFrame({'weight':weight,'return':data.loc[:i-1,'Return']}).sort_values('return')
    sorted_data.reset_index(drop=True,inplace=True)
    cum_weight = 0
    counter = 0
    while cum_weight < 0.01:
        cum_weight += sorted_data.loc[counter,'weight']
        counter += 1
    data.loc[i,'weighted_var'] = sorted_data.loc[counter - 1, 'return']

data.loc[252:,('Date','Return','history_VaR','weighted_var')].plot(x='Date')
plt.title('History VaR & Weighted VaR')
plt.show()
hist_exception = (data['Return'] < data['history_VaR'].shift(1)).value_counts()[True]
weighted_exception = (data['Return'] < data['weighted_var'].shift(1)).value_counts()[True]

def chi_test(m):
    days = len(data.loc[252:,:])
    chi_critical = chi2.ppf(0.95, 1)
    test_value = -2 * math.log(c ** (days - m) * (1 - c) ** m) + 2 * math.log((1 - m / days) ** (days - m) * (m / days) ** m)
    test_result = pd.Series({'Test Value': test_value, 'Chi-square Critical': chi_critical})
    print(test_result)
chi_test(hist_exception)
chi_test(weighted_exception)

# 2
data.loc[252:,'mean'] = data.index[252:].map(lambda x: data.loc[:x-1,'Return'].mean())
cdf_inverse = norm.ppf(0.01)
cdf = norm.ppf(0.975)
data.loc[252:,'std'] = data.index[252:].map(lambda x: data.loc[:x-1,'Return'].std())
data['f_x'] = 1 / (data['std'] * math.sqrt(2 * math.pi)) * np.exp(-(data['std'] * cdf_inverse)**2 / (2 * data['std']**2))
data.loc[252:,'cinv_upper'] = data.index[252:].map(lambda x: 1 / data.loc[x,'f_x'] * math.sqrt(0.99 * 0.01 / x) * cdf + data.loc[x,'history_VaR'])
data.loc[252:,'cinv_lower'] = data.index[252:].map(lambda x: 1 / data.loc[x,'f_x'] * math.sqrt(0.99 * 0.01 / x) * -cdf + data.loc[x,'history_VaR'])

data.loc[252:,('Date','Return','history_VaR')].plot(x='Date')
plt.fill_between(data.loc[252:,'Date'],data.loc[252:,'cinv_lower'],data.loc[252:,'cinv_upper'],color='yellow')
plt.title('Parametric interval')
plt.show()

var = np.zeros(1000)
for i in range(252,length):
    for j in range(1000):
        returns = np.random.choice(data.loc[:i-1,'Return'],size=i,replace=True)
        var[j] = np.sort(returns)[math.ceil(i * 0.01) - 1]
    data.loc[i,'bs_hist_upper'] = np.sort(var)[974]
    data.loc[i,'bs_hist_lower'] = np.sort(var)[24]


data.loc[252:,('Date','Return','history_VaR')].plot(x='Date')
#data['history_VaR'].plot()
plt.fill_between(data.loc[252:,'Date'],data.loc[252:,'bs_hist_lower'],data.loc[252:,'bs_hist_upper'],
        color='yellow')
plt.title('Bootstrap interval for historical VaR')
plt.show()

for i in range(252, length):
    weight = [lamb ** (i - j) * (1 - lamb) / (1 - lamb ** i) for j in range(1, i + 1)]
    for j in range(1000):
        indexes = np.sort(np.random.choice(data.index[:i],replace=True,size=i))
        sorted_data = pd.DataFrame({'weight':weight,'return':data.loc[indexes,'Return']}).sort_values('return')
        sorted_data.reset_index(drop=True, inplace=True)
        cum_weight = 0
        counter = 0
        while cum_weight < 0.01:
            cum_weight += sorted_data.loc[counter, 'weight']
            counter += 1
        var[j] = sorted_data.loc[counter - 1, 'return']
    data.loc[i, 'bs_wt_upper'] = np.sort(var)[974]
    data.loc[i, 'bs_wt_lower'] = np.sort(var)[24]

data.loc[252:,('Date','Return','weighted_var')].plot(x='Date')
plt.fill_between(data.loc[252:,'Date'],data.loc[252:,'bs_wt_lower'],data.loc[252:,'bs_wt_upper'],
    color='yellow')
plt.title('Bootstrap interval for weighted VaR')
plt.show()

# 3
data['monthly_std'] = data['Return'].rolling(22).std().shift(1)
data['normalized_gain'] = (data['Return'] - data['Return'].rolling(22).mean().shift(1)) / data['monthly_std']

sns.distplot(data['Return'])
plt.title('Original return distribution')
plt.show()
sns.distplot(data.loc[22:,'normalized_gain'])
plt.title('Normalized return distribution')
plt.show()

# 4
for i in range(252,length):
    weight = [lamb ** (i - j) * (1 - lamb) / (1 - lamb ** i) for j in range(1,i+1)]
    sorted_data = pd.DataFrame({'weight':weight,'return':data.loc[:i-1,'normalized_gain']}).sort_values('return')
    sorted_data.reset_index(drop=True,inplace=True)
    cum_weight = 0
    counter = 0
    while cum_weight < 0.01:
        cum_weight += sorted_data.loc[counter,'weight']
        counter += 1
    data.loc[i,'normalized_wt_var'] = sorted_data.loc[counter - 1, 'return']

normalized_exception = (data['normalized_gain'] < data['normalized_wt_var'].shift(1)).value_counts()[True]
original_skewness = data['Return'].skew()
original_kurt = data['Return'].kurt()
normalized_skewness = data['normalized_gain'].skew()
normalized_kurt = data['normalized_gain'].kurt()
print('Original skewness is',original_skewness)
print('Original kurtosis is',original_kurt)
print('Normalized skewness is',normalized_skewness)
print('Normalized kurtosis is',normalized_kurt)
print('The NO. of exceptions for normalized VaR is',normalized_exception)
