import pandas as pd
import numpy as np
from scipy.stats import skew
import matplotlib.pyplot as plt

# 1
def PS1_Q1(data):
    data['date'] = pd.to_datetime(data['date'].astype(str))
    data.loc[data['DLRET'].apply(lambda x: str(x).isalpha()), 'DLRET'] = np.nan
    data.loc[data['RET'].apply(lambda x: str(x).isalpha()), 'RET'] = np.nan
    na_return = data[data['RET'].isna()].index
    not_na_delist = data[data['DLRET'].notna()].index
    is_na = np.intersect1d(na_return, not_na_delist)
    not_na = data[data['RET'].notna() & data['DLRET'].notna()].index
    data['cum_ret'] = data['RET']
    data.loc[is_na, 'cum_ret'] = data.loc[is_na, 'DLRET']
    data.loc[not_na, 'cum_ret'] = (data.loc[not_na, 'RET'].astype(float) + 1) * (data.loc[not_na, 'DLRET'].astype(float) + 1) - 1
    data['RET'] = data['RET'].astype(float)
    data['cum_ret'] = data['cum_ret'].astype(float)
    output = data.copy()
    return output

stock = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/hw3/stocks.csv')
stock = stock.loc[((stock['EXCHCD'] == 1) | (stock['EXCHCD'] == 2) | (stock['EXCHCD'] == 3)) & ((stock['SHRCD'] == 10) | (stock['SHRCD'] == 11)), :]
stock.reset_index(inplace=True, drop=True)
stock_output = PS1_Q1(stock)
#stock['date'] = pd.to_datetime(stock['date'],format='%Y%m%d')
output1 = stock_output[['PERMNO','date','EXCHCD','cum_ret']].copy()
output1['Year'] = output1['date'].dt.year
output1['Month'] = output1.date.dt.month
output1['lag_Mkt_Cap'] = abs(stock_output['PRC']) * stock_output['SHROUT']
output1['lag_Mkt_Cap'] = output1[['PERMNO','lag_Mkt_Cap']].groupby('PERMNO').shift(1)
output1['log_return'] = np.log(1 + stock_output['RET'])
output1['Ranking_Ret'] = output1[['PERMNO','log_return']].groupby('PERMNO').rolling(11).sum()['log_return'].values
output1['Ranking_Ret'] = output1[['PERMNO','Ranking_Ret']].groupby('PERMNO').shift(2)
output1 = output1[['Year','Month','PERMNO','EXCHCD','lag_Mkt_Cap','cum_ret','Ranking_Ret']].rename(columns={'cum_ret':'RET'})
#breakpoints = output1[['Year','Month','Ranking_Ret']].groupby(['Year','Month']).apply(lambda x:pd.qcut(x['Ranking_Ret'],10,labels=False,duplicates='drop'))
print(output1)
# 2

output1_input = output1.copy()
output2 = output1_input[['Year','Month','PERMNO','lag_Mkt_Cap','RET']].copy()
output2['DM_decile'] = output1_input[['Year','Month','Ranking_Ret']].groupby(['Year','Month'])['Ranking_Ret'].transform(lambda x: pd.qcut(x,10,labels=False,duplicates='drop'))
output2['DM_decile'] = output2['DM_decile'] + 1
#df['C'] = df.groupby(['A'])['B'].transform(lambda x: pd.qcut(x, 3, labels=range(1,4)))
temp = stock_output[['date','EXCHCD','PERMNO']].copy()
temp['Ranking_Ret'] = output1_input['Ranking_Ret']
temp = temp.loc[temp['EXCHCD'] == 1,:]
breakpoints = pd.DataFrame(temp[['date','Ranking_Ret']].groupby('date')['Ranking_Ret'].apply(lambda x:pd.qcut(x,10,labels=False,duplicates='drop',retbins=True)[1]))
#breakpoints = pd.DataFrame(breakpoints)
breakpoints.reset_index(inplace=True)
breakpoints.rename(columns={'Ranking_Ret':'breakpoints'},inplace=True)
output1_input['date'] = stock_output['date']
output1_input = pd.merge(output1_input,breakpoints,on='date')
output1_input['KRF_decile'] = output1_input[['Ranking_Ret','breakpoints']].apply(lambda x: np.digitize(x[0],x[1]),axis=1)
output1_input['KRF_decile'] = output1_input['KRF_decile'].replace({11:10,0:1})

output2 = pd.merge(output2,output1_input[['Year','Month','PERMNO','KRF_decile']],on=['Year','Month','PERMNO'])
#output2['RET'] = output2['RET'].astype('float')
print(output2)
# 3
output2_input = output2.copy()

FF_mkt = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/hw3/F-F_Research_Data_Factors.CSV',skiprows=4, names=['date', 'Mkt-RF', 'SMB', 'HML', 'RF'])
FF_mkt = FF_mkt[:1110]
FF_mkt['date'] = pd.to_datetime(FF_mkt['date'],format='%Y%m')
FF_mkt['Year'] = FF_mkt['date'].dt.year
FF_mkt['Month'] = FF_mkt['date'].dt.month
years = list()
month = list()
decile = list()
for i in range(1927,2019):
    for j in range(1,13):
        for k in range(1,11):
            years.append(i)
            month.append(j)
            decile.append(k)
output3 = pd.DataFrame({'Year':years,'Month':month,'decile':decile})
#output3.apply(lambda x:)
output2_input = pd.merge(output2_input,output2_input[['Year','Month','DM_decile','lag_Mkt_Cap']].groupby(['Year','Month','DM_decile']).sum(),on=['Year','Month','DM_decile']).rename(columns={'lag_Mkt_Cap_x':'lag_Mkt_Cap','lag_Mkt_Cap_y':'DM_total_mkt'})
output2_input = pd.merge(output2_input,output2_input[['Year','Month','KRF_decile','lag_Mkt_Cap']].groupby(['Year','Month','KRF_decile']).sum(),on=['Year','Month','KRF_decile']).rename(columns={'lag_Mkt_Cap_x':'lag_Mkt_Cap','lag_Mkt_Cap_y':'KRF_total_mkt'})
DM_Ret = output2_input.groupby(['Year','Month','DM_decile'])[['RET','lag_Mkt_Cap','DM_total_mkt']].apply(lambda x:(x['lag_Mkt_Cap']/x['DM_total_mkt']*x['RET']).sum())
KRF_Ret = output2_input.groupby(['Year','Month','KRF_decile'])[['RET','lag_Mkt_Cap','KRF_total_mkt']].apply(lambda x:(x['lag_Mkt_Cap']/x['KRF_total_mkt']*x['RET']).sum())
output3['DM_Ret'] = DM_Ret.values * 100
output3['KRF_Ret'] = KRF_Ret.values * 100
output3 = pd.merge(FF_mkt[['Year','Month','RF']],output3,on=['Year','Month'])
print(output3)
# 4
output3_input = output3.copy()
output3_input['r-rf'] = output3_input['DM_Ret'] - output3_input['RF'].astype('float')
#pd.pivot_table(output3_input,values='r-rf',index=['r-rf','sigma'],columns=['decile'],aggfunc=['np.mean','np.std'])
excess_return = output3_input[['decile','r-rf']].groupby('decile').mean() * 12
standard_deviation = output3_input[['decile','r-rf']].groupby('decile').std() * np.sqrt(12)
#output4 = pd.DataFrame({'decile':[x for x in range(1,11)],'r-rf':excess_return,'sigma':standard_deviation})
Sharpe_Ratio = excess_return / standard_deviation
skewness = output3_input.groupby('decile')['DM_Ret'].apply(lambda x:np.log(x/100 + 1).skew())
wml = output3_input.loc[output3_input['decile'] == 10,'DM_Ret'].values - output3_input.loc[output3_input['decile'] == 1,'DM_Ret'].values
a = np.array([excess_return["r-rf"],standard_deviation['r-rf'],Sharpe_Ratio['r-rf'],skewness])
output4 = pd.DataFrame(a,index=['r-rf','sigma','SR','sk(m)'],columns=['Decile 1','Decile 2','Decile 3','Decile 4','Decile 5','Decile 6','Decile 7','Decile 8','Decile 9','Decile 10'])
output4['WML'] = [wml.mean() * 12, wml.std() * np.sqrt(12), wml.mean() * 12 / (wml.std() * np.sqrt(12)), skew(np.log(1 + wml / 100))]
print(output4)
# 5

DM_Ret_df = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/hw3/m_m_pt_tot.csv',header=None)
KRF_Ret_df = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/hw3/10_Portfolios_Prior_12_2.CSV',skiprows=10)[:1104]
KRF_Ret_df = KRF_Ret_df.astype(float)
DM_cor = np.zeros(10)
KRF_cor = np.zeros(10)
for i in range(1,11):
    KRF_cor[i - 1] = np.corrcoef(KRF_Ret_df.iloc[:,i],output3.loc[output3['decile'] == i,'KRF_Ret'])[0][1]
    DM_cor[i - 1] = np.corrcoef(DM_Ret_df.loc[DM_Ret_df[1] == i,2] * 100, output3.loc[(output3['decile'] == i) & (output3['Year'] < 2017),'DM_Ret'])[0][1]
output5 = pd.DataFrame(np.array([DM_cor,KRF_cor]),index=['DM correlation','KRF correlation'],columns=['Decile 1','Decile 2','Decile 3','Decile 4','Decile 5','Decile 6','Decile 7','Decile 8','Decile 9','Decile 10'])
wml_dm_cor = np.corrcoef(DM_Ret_df.loc[DM_Ret_df[1] == 10, 2].values - DM_Ret_df.loc[DM_Ret_df[1] == 1, 2].values, wml[:len(DM_Ret_df.loc[DM_Ret_df[1] == 10, 2])]/100)[0][1]
wml_krf_cor = np.corrcoef(KRF_Ret_df['Hi PRIOR'] - KRF_Ret_df['Lo PRIOR'], output3_input.loc[output3_input['decile'] == 10,'KRF_Ret'].values - output3_input.loc[output3_input['decile'] == 1,'KRF_Ret'].values)[0][1]
output5['WML'] = [wml_dm_cor,wml_krf_cor]
print(output5)

# 6

plt.plot(FF_mkt.loc[FF_mkt['date'] > '2005-12-31','date'],np.cumsum(np.log(FF_mkt.loc[FF_mkt['date'] > '2005-12-31','Mkt-RF'] / 100 + 1)),label='Market Excess Return')
plt.ylabel('Cumulative Return')
plt.plot(FF_mkt.loc[FF_mkt['date'] > '2005-12-31','date'],np.cumsum(np.log(wml[-156:]/100 + 1)),label='Momentum WML')
plt.legend()
plt.show()