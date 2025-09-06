import pandas as pd
import numpy as np
from scipy import stats
import math
#import QAM_hw1
crsp_bonds = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/hw2/crsp_bonds.csv')
crsp_bonds['MCALDT'] = pd.to_datetime(crsp_bonds['MCALDT'])
crsp_bonds['capital'] = crsp_bonds[['KYCRSPID','TMTOTOUT']].groupby('KYCRSPID').shift(1)
crsp_bonds.dropna(subset=['capital'],inplace=True)
output = crsp_bonds[['MCALDT','capital']].groupby('MCALDT').sum()
crsp_bonds.replace(-99,np.nan,inplace=True)
output['Year'] = output.index.year
output['Month'] = output.index.month
output['Bond_Ew_Ret'] = crsp_bonds[['MCALDT','TMRETNUA']].groupby('MCALDT').mean()
output.rename(columns={'capital':'Bond_lag_MV'},inplace=True)
output.reset_index(inplace=True)
crsp_bonds = pd.merge(crsp_bonds,output[['MCALDT','Bond_lag_MV']],how='left',on='MCALDT')
crsp_bonds['weight'] = crsp_bonds['capital'] / crsp_bonds['Bond_lag_MV']
output['Bond_VW_Ret'] = pd.DataFrame({'date':crsp_bonds['MCALDT'],'x':crsp_bonds['weight'] * crsp_bonds['TMRETNUA']}).groupby('date').sum().values
output.dropna(subset=['Bond_VW_Ret'],inplace=True)

#2
def PS1_Q1(data):
    data['date'] = pd.to_datetime(data['date'].astype(str))
    data['market_cap'] = np.multiply(abs(data['PRC']), data['SHROUT'])
    data['market_cap'] = data[['PERMNO', 'market_cap']].groupby('PERMNO').shift(1)
    data.dropna(subset=['market_cap'], inplace=True)
    output = data[['market_cap', 'date']].groupby('date').sum()
    data.loc[data['DLRET'].apply(lambda x: str(x).isalpha()), 'DLRET'] = np.nan
    data.loc[data['RET'].apply(lambda x: str(x).isalpha()), 'RET'] = np.nan
    na_return = data[data['RET'].isna()].index
    not_na_delist = data[data['DLRET'].notna()].index
    is_na = np.intersect1d(na_return, not_na_delist)
    not_na = data[data['RET'].notna() & data['DLRET'].notna()].index
    data.loc[is_na, 'RET'] = data.loc[is_na, 'DLRET']
    data.loc[not_na, 'RET'] = (data.loc[not_na, 'RET'].astype(float) + 1) * (data.loc[not_na, 'DLRET'].astype(float) + 1) - 1
    data.dropna(subset=['PRC', 'RET', ], inplace=True)
    data['RET'] = data['RET'].astype(float)
    output['Year'] = output.index.year
    output['Month'] = output.index.month
    output['Stock Ew Ret'] = data[['date', 'RET']].groupby('date').mean()
    output.dropna(subset=['market_cap'],inplace=True)
    output.reset_index(inplace=True)
    data = pd.merge(data, output[['date', 'market_cap']], how='left', on=['date'])
    output['Stock Vw Ret'] = pd.DataFrame({'date': data['date'], 'x': data['market_cap_x'] / data['market_cap_y'] * data['RET']}).groupby('date').sum().values
    output = output[['Year', 'Month', 'market_cap', 'Stock Ew Ret', 'Stock Vw Ret']]
    return output

stock_data = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/hw2/stocks.csv')
stock_data = stock_data.loc[((stock_data['EXCHCD'] == 1) | (stock_data['EXCHCD'] == 2) | (stock_data['EXCHCD'] == 3)) & ((stock_data['SHRCD'] == 10) | (stock_data['SHRCD'] == 11)),:]
stock_data.reset_index(inplace=True,drop=True)
output_stock = PS1_Q1(stock_data)

monthly_crsp_riskless = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/hw2/riskfree.csv')
monthly_crsp_riskless['MCALDT'] = pd.to_datetime(monthly_crsp_riskless['caldt'],format='%Y%m%d')
output_q2 = pd.merge(output,output_stock,how = 'left',on=['Year','Month'])
output_q2 = pd.merge(output_q2,monthly_crsp_riskless[['MCALDT','t90ret','t30ret']],how='left',on='MCALDT')
output_q2['Stock_Excess_VW_Ret'] = output_q2['Stock Vw Ret'] - output_q2['t30ret']
output_q2['Bond_Excess_VW_Ret'] = output_q2['Bond_VW_Ret'] - output_q2['t30ret']
output_q2 = output_q2[['Year','Month','market_cap','Stock_Excess_VW_Ret','Bond_lag_MV','Bond_Excess_VW_Ret']]
output_q2.rename(columns={'market_cap':'Stock_lag_MV'},inplace=True)
output_q2['Stock_lag_MV'] = output_q2['Stock_lag_MV'] / 1000


output_q3 = output_q2[['Year','Month','Stock_Excess_VW_Ret','Bond_Excess_VW_Ret']].copy()
output_q3['total_market_cap'] = output_q2['Bond_lag_MV'] + output_q2['Stock_lag_MV']
output_q3['Excess_Vw_Ret'] = output_q2['Bond_lag_MV'] / output_q3['total_market_cap'] * output_q3['Bond_Excess_VW_Ret'] + output_q2['Stock_lag_MV'] / output_q3['total_market_cap'] * output_q3['Stock_Excess_VW_Ret']
output_q3['Excess_60_40_Ret'] = 0.6 * output_q3['Stock_Excess_VW_Ret'] + 0.4 * output_q3['Bond_Excess_VW_Ret']
output_q3['Stock_inverse_sigma_hat'] = 1 / output_q3['Stock_Excess_VW_Ret'].rolling(36).std().shift(1)
output_q3['Bond_inverse_sigma_hat'] = 1 / output_q3['Bond_Excess_VW_Ret'].rolling(36).std().shift(1)
output_q3['Unlevered_k'] = 1 / (output_q3['Stock_inverse_sigma_hat'] + output_q3['Bond_inverse_sigma_hat'])
output_q3['Excess_Unlevered_RP_Ret'] = (output_q3['Unlevered_k'] * output_q3['Stock_inverse_sigma_hat']).shift(1) * output_q3['Stock_Excess_VW_Ret'] + (output_q3['Unlevered_k'] * output_q3['Bond_inverse_sigma_hat']).shift(1) * output_q3['Bond_Excess_VW_Ret']
constant_k = output_q3.loc[output_q3['Year'] > 1928,'Excess_Vw_Ret'].std() / ((output_q3['Stock_inverse_sigma_hat']).shift(1) * output_q3['Stock_Excess_VW_Ret'] + (output_q3['Bond_inverse_sigma_hat']).shift(1) * output_q3['Bond_Excess_VW_Ret']).std()
output_q3['Excess_Levered_RP_Ret'] = (constant_k * output_q3['Stock_inverse_sigma_hat']).shift(1) * output_q3['Stock_Excess_VW_Ret'] + (constant_k * output_q3['Bond_inverse_sigma_hat']).shift(1) * output_q3['Bond_Excess_VW_Ret']

index = ['CRSP stocks','CRSP bonds','Value-weighted portfolio','60/40 portfolio','unlevered RP','levered RP']
q4_data = output_q3.loc[48:1014,('Stock_Excess_VW_Ret','Bond_Excess_VW_Ret','Excess_Vw_Ret','Excess_60_40_Ret','Excess_Unlevered_RP_Ret','Excess_Levered_RP_Ret')]
observations = len(q4_data)
Annualized_Mean = q4_data.mean() * 12
t_stats = q4_data.apply(lambda x: x.mean()/x.std() * math.sqrt(observations))
Annualized_std = q4_data.std() * math.sqrt(12)
Annualized_SR = Annualized_Mean / Annualized_std
Skewness = q4_data.skew()
Excess_Kurtosis = q4_data.kurt()
output_q4 = pd.DataFrame({'Annualized Mean':Annualized_Mean.values,'T stat':t_stats.values,'Annualized Standard Deviation':Annualized_std.values,'Annualized Sharpe Ratio':Annualized_SR.values,'Skewness':Skewness.values,'Excess Kurtosis':Excess_Kurtosis.values},index=index)