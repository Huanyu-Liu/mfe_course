#%%
import pandas as pd
import numpy as np
import wrds
import scipy.stats as sts
import statsmodels.api as sm
import scipy.stats as ss
from pandas.tseries.offsets import *
import matplotlib.pyplot as plt
#%%
q_factors = pd.read_excel('/Users/huanyu/Desktop/QuantitativeAssetManagement/final/HXZ q-Factors (monthly 1967 to 2018).xlsx')
q_factors = q_factors.loc[(q_factors['Year'] > 1971) & (q_factors['Year'] < 2019),:].reset_index(drop=True)
q_factors[['RF','MKT-RF','ME','IA','ROE']] = q_factors[['RF','MKT-RF','ME','IA','ROE']] / 100
q_factors['MKT'] = q_factors['MKT-RF'] + q_factors['RF']
ff_3factor = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/final/F-F_Research_Data_Factors_CSV.zip',skiprows=3).rename(columns={'Unnamed: 0':'date'})[:1113].astype(float)
ff_3factor['date'] = pd.to_datetime(ff_3factor['date'],format='%Y%m')
ff_3factor = ff_3factor.loc[(ff_3factor['date'] < '2019-01-01') & (ff_3factor['date'] > '1971-12-01'),:].reset_index(drop=True)
ff_3factor[['Mkt-RF','SMB','HML','RF']] = ff_3factor[['Mkt-RF','SMB','HML','RF']] / 100
ff_mom = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/final/F-F_Momentum_Factor_CSV.zip',skiprows=14,names=['date','mom'])[:1107].astype(float)
ff_mom['date'] = pd.to_datetime(ff_mom['date'],format='%Y%m')
ff_mom = ff_mom.loc[(ff_mom['date'] < '2019-01-01') & (ff_mom['date'] > '1971-12-01'),:].reset_index(drop=True)
ff_mom['mom'] /= 100
#%%
means = q_factors[['ME','IA','ROE']].mean()
sts.ttest_1samp(q_factors[['ME','IA','ROE']],0)
# panel A
def stats_func(y,x):
    lm = sm.OLS(y,sm.add_constant(x))
    res = lm.fit()
    return res.params, res.tvalues, res.rsquared
stats_func(q_factors['ME'],sm.add_constant(q_factors['MKT-RF']))
stats_func(q_factors['ME'],sm.add_constant(ff_3factor[['Mkt-RF','SMB','HML']]))
stats_func(q_factors['ME'],sm.add_constant(pd.concat([ff_3factor[['Mkt-RF','SMB','HML']],ff_mom['mom']], axis=1)))
stats_func(q_factors['IA'],sm.add_constant(q_factors['MKT-RF']))
stats_func(q_factors['IA'],sm.add_constant(ff_3factor[['Mkt-RF','SMB','HML']]))
stats_func(q_factors['IA'],sm.add_constant(pd.concat([ff_3factor[['Mkt-RF','SMB','HML']],ff_mom['mom']], axis=1)))
stats_func(q_factors['ROE'],sm.add_constant(q_factors['MKT-RF']))
stats_func(q_factors['ROE'],sm.add_constant(ff_3factor[['Mkt-RF','SMB','HML']]))
stats_func(q_factors['ROE'],sm.add_constant(pd.concat([ff_3factor[['Mkt-RF','SMB','HML']],ff_mom['mom']], axis=1)))

# panel B
question1_df = pd.concat([q_factors[['ME','IA','ROE','MKT-RF']],ff_3factor[['SMB','HML']],ff_mom['mom']], axis=1)
n = len(question1_df)
panel_B = question1_df.corr()
t = panel_B * np.sqrt((n-2)/(1 - panel_B * panel_B))
p_values = ss.t.sf(np.abs(t), n-1) * 2
t=panel_B / np.sqrt((1-panel_B*panel_B) / (n - 2))

#%%
conn = wrds.Connection()
comp_annual = conn.raw_sql("""
                    select gvkey, at, datadate, fyear
                    from compa.funda
                    where indfmt='INDL'
                    and datafmt='STD'
                    and popsrc='D'
                    and consol='C'
                    and datadate >= '01/01/1970'
                    """)
#comp_annual.to_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/final/comp_annual.csv',index=False)
comp_quarter = conn.raw_sql("""
                    select gvkey, datadate, fqtr, fyearq, atq, txditcq,
                    seqq, ceqq, pstkq, ltq, pstkrq, ibq, rdq, cstkcvq
                    from compa.fundq
                    where indfmt='INDL'
                    and datafmt='STD'
                    and popsrc='D'
                    and consol='C'
                    and datadate >= '01/01/1965'
                    """)
# comp_quarter.to_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/final/comp_quarter.csv',index=False)
#%%
# comp_annual = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/final/comp_annual.csv')
# comp_quarter = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/final/comp_quarter.csv')
comp_annual['gvkey'] = comp_annual['gvkey'].astype(int)
comp_annual['datadate'] = pd.to_datetime(comp_annual['datadate'])
comp_annual['year'] = comp_annual['datadate'].dt.year
comp_annual = comp_annual.sort_values(['gvkey','datadate'])
comp_annual['lag_at'] = comp_annual[['gvkey','at']].groupby('gvkey').shift(1)
comp_annual['IA'] = comp_annual['at'] / comp_annual['lag_at'] - 1
#%%
comp_quarter['gvkey'] = comp_quarter['gvkey'].astype(int)
comp_quarter['datadate'] = pd.to_datetime(comp_quarter['datadate'])
comp_quarter['year'] = comp_quarter['datadate'].dt.year
comp_quarter['rdq'] = pd.to_datetime(comp_quarter['rdq'])
comp_quarter['she'] = np.where(comp_quarter['seqq'].isnull(),comp_quarter['ceqq'] + comp_quarter['pstkq'],comp_quarter['seqq'])
comp_quarter['she'] = np.where(comp_quarter['seqq'].isnull(),comp_quarter['atq'] - comp_quarter['ltq'],comp_quarter['she'])
comp_quarter['she'] = comp_quarter['she'].fillna(0)
comp_quarter['txditcq'] = comp_quarter['txditcq'].fillna(0)
comp_quarter['ps'] = np.where(comp_quarter['pstkrq'].isnull(), comp_quarter['cstkcvq'], comp_quarter['pstkrq'])
comp_quarter['ps'] = comp_quarter['ps'].fillna(0)
comp_quarter['be'] = comp_quarter['she'] + comp_quarter['txditcq'] - comp_quarter['ps']
comp_quarter['be'] = np.where(comp_quarter['be'] > 0, comp_quarter['be'], np.nan)
comp_quarter = comp_quarter.sort_values(['gvkey','datadate'])

comp_quarter['lag_be'] = comp_quarter[['gvkey','be']].groupby('gvkey').shift(1)

comp_quarter['ROE'] = comp_quarter['ibq'] / comp_quarter['lag_be']
comp_quarter = comp_quarter.dropna(subset=['be','lag_be']).reset_index(drop=True)
comp_quarter['jdate'] = comp_quarter['rdq'] + MonthEnd(2)

#%%
crsp_m = conn.raw_sql("""
                      select a.permno, a.permco, a.date, b.shrcd, b.exchcd,
                      a.ret, a.shrout, a.prc
                      from crspa.msf as a
                      left join crspa.msenames as b
                      on a.permno=b.permno
                      and b.namedt<=a.date
                      and a.date<=b.nameendt
                      where a.date between '01/01/1970' and '12/31/2018'
                      and b.exchcd between 1 and 3
                      and b.shrcd between 10 and 11
                      """)
# crsp_m.to_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/final/crsp_m.csv',index=False)
#%%
# crsp_m = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/final/crsp_m.csv')
crsp_m[['permco','permno','shrcd','exchcd']]=crsp_m[['permco','permno','shrcd','exchcd']].astype(int)
crsp_m['date'] = pd.to_datetime(crsp_m['date'])
crsp_m['jdate'] = crsp_m['date'] + MonthEnd(0) ## month end 0 month end 1
crsp_m['ret']= crsp_m['ret'].fillna(0)
#%%
dlret = conn.raw_sql("""
                     select permno, dlret, dlstdt 
                     from crsp.msedelist
                     """)
# dlret.to_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/final/dlret.csv',index=False)
#%%
# dlret = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/final/dlret.csv')
dlret['dlstdt'] = pd.to_datetime(dlret['dlstdt'])
dlret['jdate'] = dlret['dlstdt'] + MonthEnd(0)
dlret['permno'] = dlret['permno'].astype(int)
# dlret['dlret'] = dlret['dlret'].fillna(0)
crsp = pd.merge(crsp_m, dlret, how='left', on=['permno','jdate'])
crsp['dlret'] = crsp['dlret'].fillna(0)
crsp['adj_ret'] = (1 + crsp['dlret']) * (1 + crsp['ret']) - 1
crsp['me'] = np.abs(crsp['prc']) * crsp['shrout']
crsp = crsp.sort_values(['jdate','permco'])
max_me = crsp[['jdate','permco','me']].groupby(['jdate','permco']).max().reset_index().rename(columns={'me':'max_me'})
crsp1 = pd.merge(crsp,max_me,how='left',on=['jdate','permco'])
crsp1 = crsp1.sort_values(['permno','jdate'])
crsp1['ffdate'] = crsp1['jdate'] + MonthEnd(-6)
crsp1['ffyear'] = crsp1['ffdate'].dt.year
crsp1['ffmonth'] = crsp1['ffdate'].dt.month
crsp1['lag_me'] = crsp1[['permno','max_me']].groupby('permno').shift(1)
#%%
ccm=conn.raw_sql("""
                  select gvkey, lpermno as permno, linktype, linkprim,
                  linkdt, linkenddt
                  from crsp.ccmxpf_linktable
                  where substr(linktype,1,1)='L'
                  and (linkprim ='C' or linkprim='P')
                  """)
# ccm.to_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/final/ccm.csv',index=False)
#%%
# ccm = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/final/ccm.csv')
ccm[['gvkey','permno']] = ccm[['gvkey','permno']].astype(int)
ccm['linkdt'] = pd.to_datetime(ccm['linkdt'])
ccm['linkenddt'] = pd.to_datetime(ccm['linkenddt'])
ccm['linkenddt'] = ccm['linkenddt'].fillna(pd.to_datetime('today'))
ccm1 = pd.merge(comp_annual[['gvkey','datadate','IA']],ccm,how='left',on=['gvkey'])
ccm1['jdate'] = ccm1['datadate'] + YearEnd(0) + MonthEnd(6)
ccm2 = ccm1.loc[(ccm1['jdate'] >= ccm1['linkdt']) & (ccm1['jdate'] <= ccm1['linkenddt']),:].copy()
ccm2['ffyear'] = ccm2['jdate'].dt.year
ccm_IA = pd.merge(crsp1,ccm2[['gvkey','permno','ffyear','IA']],how='left',on=['permno','ffyear'])
ccm3 = pd.merge(comp_quarter[['gvkey','jdate','ROE']],ccm,how='left',on=['gvkey'])
ccm4 = ccm3.loc[(ccm3['jdate'] >= ccm3['linkdt']) & (ccm3['jdate'] <= ccm3['linkenddt']),:].copy()
ccm_ROE = pd.merge(crsp1[['permno','jdate']],ccm4,how='left',on=['permno','jdate'])
ccm_ROE['ROE'] = ccm_ROE['ROE'].fillna(method='ffill')
ccm_final = pd.merge(ccm_IA,ccm_ROE[['permno','jdate','ROE']],how='left',on=['permno','jdate'])
ccm_final = ccm_final.loc[(ccm_final['jdate'] > '1972-12-31') & (ccm_final['jdate'] < '2019-01-01'),:].dropna(subset=['IA','ROE','lag_me']).reset_index(drop=True)
#%%
ccm_july = ccm_final.loc[ccm_final['ffmonth'] == 1,:].reset_index(drop=True)
ccm_NYSE = ccm_july.loc[ccm_july['exchcd'] == 1].reset_index(drop=True)
size_bins = ccm_NYSE.groupby('ffyear')['lag_me'].apply(lambda x: pd.qcut(x,2,retbins=True)[1][1:-1])
#size_bins
ccm_july['size_decile'] = ccm_july[['ffyear','lag_me']].apply(lambda x: np.digitize(x['lag_me'],size_bins[x['ffyear']]),axis=1)
#ccm_july['decile'] = ccm_july.groupby('ffyear')['lag_me'].apply(lambda x: x.name)
IA_bins = ccm_NYSE.groupby('ffyear')['IA'].apply(lambda x: pd.qcut(x, [0,0.3,0.7,1], retbins=True)[1][1:-1])
ccm_july['IA_decile'] = ccm_july[['ffyear','IA']].apply(lambda x: np.digitize(x['IA'],IA_bins[x['ffyear']]),axis=1)
ccm_NYSE_month = ccm_final.loc[(ccm_final['exchcd'] == 1) & (ccm_final['ffyear'] > 1973),:].reset_index(drop=True)
ccm_NYSE_month['year_month'] = ccm_NYSE_month['ffyear'] * 100 + ccm_NYSE_month['ffmonth']
ROE_bins = ccm_NYSE_month.groupby('year_month')['ROE'].apply(lambda x: pd.qcut(x, [0,0.3,0.7,1], retbins=True)[1][1:-1])
ROE_decile = ccm_final.loc[ccm_final['ffyear'] > 1972,['ffyear','ffmonth','ROE','permno']].reset_index(drop=True)
ROE_decile['year_month'] = ROE_decile['ffyear'] * 100 + ROE_decile['ffmonth']
ROE_decile['ROE_decile'] = ROE_decile[['year_month','ROE']].apply(lambda x: np.digitize(x['ROE'],ROE_bins[x['year_month']]),axis=1)
#%%
ccm_final = ccm_final.loc[ccm_final['ffyear'] > 1972,:].reset_index(drop=True)
ccm_decile = pd.merge(ccm_final[['ffyear','permno','gvkey','ffmonth','lag_me','adj_ret']],ccm_july[['ffyear','permno','size_decile','IA_decile']],how='left',on=['ffyear','permno'])
ccm_decile1 = pd.merge(ccm_decile,ROE_decile[['ffyear','ffmonth','permno','ROE_decile']],how='left',on=['ffyear','ffmonth','permno'])

ccm_decile1 = ccm_decile1.dropna(subset={'size_decile','IA_decile','ROE_decile','lag_me'}).reset_index(drop=True)
ccm_decile1[['size_decile','IA_decile','ROE_decile']] = ccm_decile1[['size_decile','IA_decile','ROE_decile']].astype(int)
# ccm_decile1['decile'] = ccm_decile1[['size_decile','IA_decile','ROE_decile']]
ccm_decile1['decile'] = ccm_decile1['size_decile'].astype(str) + ccm_decile1['IA_decile'].astype(str) + ccm_decile1['ROE_decile'].astype(str)
#%%
value_weight = lambda x: np.average(x, weights=ccm_decile1.loc[x.index, "lag_me"])
portfolio_func = {'adj_ret':value_weight}
portfolio_ret = ccm_decile1.groupby(['ffyear','ffmonth','decile']).agg(portfolio_func).reset_index()
mkt_ret = ccm_decile1.groupby(['jdate']).agg(portfolio_func).reset_index()
mkt_ret.columns = ['jdate', 'MKT']
portfolio_ret_pivot = portfolio_ret.pivot_table('adj_ret', ['ffyear','ffmonth'], 'decile').reset_index()
portfolio_ret_pivot['ME'] = np.nanmean([portfolio_ret_pivot['000'], portfolio_ret_pivot['001'], portfolio_ret_pivot['002'], portfolio_ret_pivot['010'], portfolio_ret_pivot['011'], portfolio_ret_pivot['012'], portfolio_ret_pivot['120'], portfolio_ret_pivot['121'], portfolio_ret_pivot['022']],axis=0) - np.nanmean([portfolio_ret_pivot['100'], portfolio_ret_pivot['101'], portfolio_ret_pivot['102'], portfolio_ret_pivot['110'], portfolio_ret_pivot['111'], portfolio_ret_pivot['112'], portfolio_ret_pivot['120'], portfolio_ret_pivot['121'], portfolio_ret_pivot['122']],axis=0)
portfolio_ret_pivot['IA'] = np.nanmean([portfolio_ret_pivot['000'], portfolio_ret_pivot['001'], portfolio_ret_pivot['002'], portfolio_ret_pivot['100'], portfolio_ret_pivot['101'], portfolio_ret_pivot['102']], axis=0) - np.nanmean([portfolio_ret_pivot['020'], portfolio_ret_pivot['021'], portfolio_ret_pivot['022'], portfolio_ret_pivot['120'], portfolio_ret_pivot['121'], portfolio_ret_pivot['122']],axis=0)
portfolio_ret_pivot['ROE'] = np.nanmean([portfolio_ret_pivot['002'], portfolio_ret_pivot['012'], portfolio_ret_pivot['022'], portfolio_ret_pivot['102'], portfolio_ret_pivot['112'], portfolio_ret_pivot['122']],axis=0) - np.nanmean([portfolio_ret_pivot['000'], portfolio_ret_pivot['010'], portfolio_ret_pivot['020'], portfolio_ret_pivot['100'], portfolio_ret_pivot['110'], portfolio_ret_pivot['120']],axis=0)

replicated_factor = pd.merge(mkt_ret, portfolio_ret_pivot, how='left', on=['jdate'])
q_factors.loc[(q_factors['Year'] < 2018) & (q_factors['Year'] > 1973),'IA'].corr(portfolio_ret_pivot.loc[(portfolio_ret_pivot['ffyear'] < 2018),'IA'])

# q_factors['date'] = pd.to_datetime(q_factors['Year'] * 100 + q_factors['Month'],format='%Y%m')
# q_factors['replcated_mkt-rf'] = q_factors['MKT-RF'] + np.random.normal(0,0.0025,564)
# q_factors.loc[q_factors['date']]
# q_factors[['date','MKT-RF']].plot(x='date')
# plt.show()
# q_factors[['date','replcated_mkt-rf']].plot(x='date')
# plt.show()
# plt.plot(q_factors['date'],np.abs(q_factors['replcated_mkt-rf'] - q_factors['MKT-RF']))
# plt.show()
replicated_factor.columns = ['jdate', 'replicated_MKT', 'replicated_ME', 'replicated_IA', 'replicated_ROE']
replicated_factor[['replicated_MKT', 'replicated_ME', 'replicated_IA', 'replicated_ROE']] *= 100
q_factors.columns = ['Year', 'Month', 'Original_RF', 'Original_MKT-RF', 'Original_ME', 'Original_IA', 'Original_ROE', 'Original_MKT']
Q2_output = pd.concat([q_factors, replicated_factor], axis=1)
Q2_output['replicated_MKT-RF'] = Q2_output['replicated_MKT'] - Q2_output['Original_RF']
Q2_output = Q2_output[['jdate', 'Original_MKT-RF', 'replicated_MKT-RF', 'Original_ME', 'replicated_ME', 'Original_IA', 'replicated_IA', 'Original_ROE', 'replicated_ROE']]
Q2_panel_A = Q2_output[['Original_MKT-RF', 'replicated_MKT-RF', 'Original_ME', 'replicated_ME', 'Original_IA', 'replicated_IA', 'Original_ROE', 'replicated_ROE']].apply(lambda x:pd.Series({'Average':x.mean()*12, 'Volatility':x.std()*np.sqrt(12), 'Sharpe Ratio': np.sqrt(12)*x.mean()/x.std(), 'Skewness':x.skew(), 'Kurtosis':x.kurtosis()}))
Q2_output['Diff_MKT-RF'] = 100 * np.abs(Q2_output['replicated_MKT-RF'] - Q2_output['Original_MKT-RF'])
Q2_output['Diff_ME'] = 100 * np.abs(Q2_output['replicated_ME'] - Q2_output['Original_ME'])
Q2_output['Diff_IA'] = 100 * np.abs(Q2_output['replicated_IA'] - Q2_output['Original_IA'])
Q2_output['Diff_ROE'] = 100 * np.abs(Q2_output['replicated_ROE'] - Q2_output['Original_ROE'])
Q2_panel_B = Q2_output[['Diff_MKT-RF', 'Diff_ME', 'Diff_IA', 'Diff_ROE']].apply(lambda x:pd.Series({'Average':x.mean(), 'Volatility':x.std(), 'Min': x.min(), '25%': x.quantile(q=0.25), '50%': x.quantile(q=0.5), '75%': x.quantile(q=0.75),  'Max':x.max()}))
Q2_panel_C = pd.DataFrame(index=['correlation'],columns=['MKT-RF','ME', 'IA', 'ROE'])
Q2_panel_C['MKT-RF'] = Q2_output['replicated_MKT-RF'].corr(Q2_output['Original_MKT-RF'])
Q2_panel_C['ME'] = Q2_output['replicated_ME'].corr(Q2_output['Original_ME'])
Q2_panel_C['IA'] = Q2_output['replicated_IA'].corr(Q2_output['Original_IA'])
Q2_panel_C['ROE'] = Q2_output['replicated_ROE'].corr(Q2_output['Original_ROE'])

plt.plot(Q2_output['Original_MKT-RF'])
plt.show()
plt.plot(Q2_output['replicated_MKT-RF'])
plt.show()
plt.plot(Q2_output['Diff_MKT-RF'])
plt.show()

plt.plot(Q2_output['Original_ME'])
plt.show()
plt.plot(Q2_output['replicated_ME'])
plt.show()
plt.plot(Q2_output['Diff_ME'])
plt.show()

plt.plot(Q2_output['Original_IA'])
plt.show()
plt.plot(Q2_output['replicated_IA'])
plt.show()
plt.plot(Q2_output['Diff_IA'])
plt.show()

plt.plot(Q2_output['Original_ROE'])
plt.show()
plt.plot(Q2_output['replicated_ROE'])
plt.show()