import pandas as pd
import numpy as np
import math
# 1

def PS1_Q1(data):
    data.loc[data['DLRET'].apply(lambda x: str(x).isalpha()), 'DLRET'] = np.nan
    data.loc[data['RET'].apply(lambda x: str(x).isalpha()), 'RET'] = np.nan
    na_return = data[data['RET'].isna()].index
    not_na_delist = data[data['DLRET'].notna()].index
    is_na = np.intersect1d(na_return, not_na_delist)
    not_na = data[data['RET'].notna() & data['DLRET'].notna()].index
    data.iloc[is_na, 6] = data.iloc[is_na, 4]
    data.iloc[not_na, 6] = (data.iloc[not_na, 6].astype(float) + 1) * (data.iloc[not_na, 4].astype(float) + 1) - 1
    data.dropna(subset=['PRC', 'RET', ], inplace=True)
    data['date'] = pd.to_datetime(data['date'].astype(str))
    data['RET'] = data['RET'].astype(float)
    data['market_cap'] = np.multiply(abs(data['PRC']), data['SHROUT'])
    data['market_cap'] = data[['PERMNO', 'market_cap']].groupby('PERMNO').shift(1)
    data.dropna(subset=['market_cap'], inplace=True)
    output = data[['market_cap', 'date']].groupby('date').sum()
    output['Year'] = output.index.year
    output['Month'] = output.index.month
    output['Stock Ew Ret'] = data[['date', 'RET']].groupby('date').mean()
    output.reset_index(inplace=True)
    data = pd.merge(data, output[['date', 'market_cap']], how='left', on=['date'])
    output['Stock Vw Ret'] = pd.DataFrame(
        {'date': data['date'], 'x': data['market_cap_x'] / data['market_cap_y'] * data['RET']}).groupby(
        'date').sum().values
    output = output[['Year', 'Month', 'market_cap', 'Stock Ew Ret', 'Stock Vw Ret']]
    return output

raw_data = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/hw1/stocks.csv')
raw_data = raw_data.loc[((raw_data['EXCHCD'] == 1) | (raw_data['EXCHCD'] == 2) | (raw_data['EXCHCD'] == 3)) & ((raw_data['SHRCD'] == 10) | (raw_data['SHRCD'] == 11)),:]
raw_data.reset_index(inplace=True)
output = PS1_Q1(raw_data)


# 2

def PS1_Q2(FF_market,output):
    FF_market['Est_Ret'] = output.loc[5:, 'Stock Vw Ret'].values * 100 - FF_market['RF'].values.astype(float)
    x = lambda a: [a.mean() / 100 * 12, a.std() / 100 * math.sqrt(12), a.mean() / a.std() * math.sqrt(12), a.skew(), a.kurtosis()]
    output2 = pd.DataFrame(
        {'Est_Excess_Ret': x(FF_market['Est_Ret']), 'Act_Excess_Ret': x(FF_market['Mkt-RF'].astype(float))},
        index=['Annual_Mean', 'Annual_sd', 'Annual_SR', 'Skewness', 'Kurtosis'])
    return output2

FF_mkt = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/hw1/F-F_Research_Data_Factors.CSV', skiprows=4, names=['date','Mkt-RF','SMB','HML','RF'])
FF_mkt = FF_mkt.iloc[FF_mkt[FF_mkt['date'] == '192607'].index[0]:FF_mkt[FF_mkt['date'] == '201812'].index[0] + 1,]
output2 = PS1_Q2(FF_mkt,output)

# 3

def PS1_Q3(FF_market, output):
    FF_market['Est_Ret'] = output.loc[5:, 'Stock Vw Ret'].values * 100 - FF_market['RF'].values.astype(float)
    correlation = np.corrcoef(FF_market['Mkt-RF'].astype(float), FF_market['Est_Ret'])[1][0]
    max_abs_diff = abs(np.subtract(FF_market['Mkt-RF'].astype(float), FF_market['Est_Ret'])).max() / 100
    return [correlation, max_abs_diff]
output3 = PS1_Q3(FF_mkt, output)