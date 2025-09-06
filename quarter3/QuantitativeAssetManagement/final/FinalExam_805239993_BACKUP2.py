
import pandas as pd
import numpy as np
import wrds
from pandas.tseries.offsets import *
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import skew, kurtosis
import scipy.stats as sts
import statsmodels.api as sm


###################
# Connect to WRDS #
###################
conn = wrds.Connection()

###################
# Compustat Block #
###################

# # Download Annual Compustat Data
# comp_A = conn.raw_sql("""
#                     select gvkey, at, datadate, fyear
#                     from compa.funda
#                     where indfmt='INDL'
#                     and datafmt='STD'
#                     and popsrc='D'
#                     and consol='C'
#                     and datadate >= '01/01/1965'
#                     """)
# comp_A.to_csv('/Users/leonard/Desktop/MFE-QAM/PS/comp_A.csv')

comp_A = pd.read_csv('/Users/leonard/Desktop/MFE-QAM/PS/comp_A.csv')

comp_A['datadate']=pd.to_datetime(comp_A['datadate']) #convert datadate to date fmt
comp_A['year']=comp_A['datadate'].dt.year

comp_A=comp_A.sort_values(by=['gvkey','datadate'])

# lagged total asset
comp_A['lat'] = comp_A.groupby(['gvkey'])['at'].shift(1)
# I/A
comp_A['IA'] = (comp_A['at'] - comp_A['lat'])/comp_A['lat']


# # Download Quarterly Compustat Data
# comp_Q = conn.raw_sql("""
#                     select gvkey, datadate, fqtr, fyearq, atq, txditcq,
#                     seqq, ceqq, pstkq, ltq, pstkrq, ibq, rdq, cstkcvq
#                     from compa.fundq
#                     where indfmt='INDL'
#                     and datafmt='STD'
#                     and popsrc='D'
#                     and consol='C'
#                     and datadate >= '01/01/1965'
#                     """)

comp_Q = pd.read_csv('/Users/leonard/Desktop/MFE-QAM/PS/comp_Q.csv')

comp_Q['datadate']=pd.to_datetime(comp_Q['datadate']) #convert datadate to date fmt
comp_Q['year']=comp_Q['datadate'].dt.year

# create Shareholders' equity
comp_Q['she']=np.where(comp_Q['seqq'].isnull(), comp_Q['ceqq']+comp_Q['pstkq'], comp_Q['seqq'])
comp_Q['she']=np.where(comp_Q['she'].isnull(), comp_Q['atq']-comp_Q['ltq'], comp_Q['she'])
comp_Q['she']=np.where(comp_Q['she'].isnull(), 0, comp_Q['she'])

# create balance-sheet deferred taxes and investment tax credit (item TXDITCQ)
comp_Q['txditcq']=np.where(comp_Q['txditcq'].isnull(), 0, comp_Q['txditcq'])

# create preferrerd stock
comp_Q['ps']=np.where(comp_Q['pstkrq'].isnull(), comp_Q['cstkcvq'], comp_Q['pstkrq'])
comp_Q['ps']=np.where(comp_Q['ps'].isnull(),0,comp_Q['ps'])

# create book equity, nan if smaller than 0
comp_Q['be']=comp_Q['she']+comp_Q['txditcq']-comp_Q['ps']
comp_Q['be']=np.where(comp_Q['be']>0, comp_Q['be'], np.nan)

comp_Q=comp_Q.sort_values(by=['gvkey','datadate'])
comp_Q['lbe'] = comp_Q.groupby(['gvkey'])['be'].shift(1)
comp_Q['ROE'] = comp_Q['ibq']/comp_Q['lbe']

# delete negative book equity
comp_Q = comp_Q.loc[comp_Q['be'].notnull()]
comp_Q = comp_Q.loc[comp_Q['lbe'].notnull()]

# Create announcing date year and month, use them to merge with CRSP
comp_Q['rdq'] = pd.to_datetime(comp_Q['rdq'])
# Add one month to merge COMP t-1 month with CRSP t month
comp_Q['jdate'] = comp_Q['rdq'] + MonthEnd(0)
comp_Q['jdate'] = comp_Q['jdate'] + MonthEnd(1)

###################
# CRSP Block      #
###################
# crsp_m = conn.raw_sql("""
#                       select a.permno, a.permco, a.date, b.shrcd, b.exchcd,
#                       a.ret, a.shrout, a.prc
#                       from crspa.msf as a
#                       left join crspa.msenames as b
#                       on a.permno=b.permno
#                       and b.namedt<=a.date
#                       and a.date<=b.nameendt
#                       where a.date between '01/01/1970' and '12/31/2018'
#                       and b.exchcd between 1 and 3
#                       and b.shrcd between 10 and 11
#                       """)
#
# crsp_m.to_csv('/Users/leonard/Desktop/MFE-QAM/PS/crsp_m.csv')

crsp_m = pd.read_csv('/Users/leonard/Desktop/MFE-QAM/PS/crsp_m.csv')

# change variable format to int
crsp_m[['permco','permno','shrcd','exchcd']]=crsp_m[['permco','permno','shrcd','exchcd']].astype(int)

# Line up date to be end of month
crsp_m['date']=pd.to_datetime(crsp_m['date'])
crsp_m['jdate']=crsp_m['date']+MonthEnd(0)

# add delisting return
dlret = conn.raw_sql("""
                     select permno, dlret, dlstdt 
                     from crsp.msedelist
                     """)
dlret.permno=dlret.permno.astype(int)
dlret['dlstdt']=pd.to_datetime(dlret['dlstdt'])
dlret['jdate']=dlret['dlstdt']+MonthEnd(0)

crsp = pd.merge(crsp_m, dlret, how='left',on=['permno','jdate'])
crsp['dlret']=crsp['dlret'].fillna(0)
crsp['ret']=crsp['ret'].fillna(0)
crsp['retadj']=(1+crsp['ret'])*(1+crsp['dlret'])-1 # cumulative return
crsp['me']=crsp['prc'].abs()*crsp['shrout'] # calculate market equity
crsp=crsp.drop(['dlret','dlstdt','prc','shrout'], axis=1)
crsp=crsp.sort_values(by=['jdate','permco','me'])

### Aggregate Market Cap ###
# maximum of me across different permno belonging to same permco a given date
crsp_summe = crsp.groupby(['jdate','permco'])['me'].max().reset_index()
# merge permco me to permno me
crsp=crsp.drop(['me'], axis=1)
crsp2=pd.merge(crsp, crsp_summe, how='left', on=['jdate','permco'])
# sort by permno and date
crsp2=crsp2.drop(['permco','date','ret'], axis=1)
crsp2=crsp2.sort_values(by=['permno','jdate'])

# t July to t+1 June is within a ffyear
crsp2['ffdate']=crsp2['jdate']+MonthEnd(-6)
crsp2['ffyear']=crsp2['ffdate'].dt.year
crsp2['ffmonth']=crsp2['ffdate'].dt.month

# lag market cap
crsp2['lme']=crsp2.groupby(['permno'])['me'].shift(1)
crsp2=crsp2.drop(['me'], axis=1)


#######################
# CCM Block           #
#######################
ccm=conn.raw_sql("""
                  select gvkey, lpermno as permno, linktype, linkprim, 
                  linkdt, linkenddt
                  from crsp.ccmxpf_linktable
                  where substr(linktype,1,1)='L'
                  and (linkprim ='C' or linkprim='P')
                  """)

ccm['linkdt']=pd.to_datetime(ccm['linkdt'])
ccm['linkenddt']=pd.to_datetime(ccm['linkenddt'])
# if linkenddt is missing then set to today date
ccm['linkenddt']=ccm['linkenddt'].fillna(pd.to_datetime('today'))

ccm['gvkey'] = ccm['gvkey'].astype('int')

# link CRSP and COMP Yearly data
ccm1=pd.merge(comp_A[['gvkey','datadate','IA']],ccm,how='left',on=['gvkey'])
ccm1['yearend']=ccm1['datadate']+YearEnd(0)
ccm1['jdate']=ccm1['yearend']+MonthEnd(6)

# set link date bounds
ccm2=ccm1[(ccm1['jdate']>=ccm1['linkdt'])&(ccm1['jdate']<=ccm1['linkenddt'])].copy()
# create ffyear and delete jdate, because jdate is not actual date any more
ccm2['ffyear']=ccm2['jdate'].dt.year
ccm2=ccm2[['gvkey','permno','ffyear','IA']]

# link comp and crsp using permno and ffyear
ccm_IA=pd.merge(crsp2, ccm2, how='left', on=['permno', 'ffyear'])
ccm_IA=ccm_IA.drop(['gvkey'],axis=1)


# link CRSP and COMP Quarterly data
ccm3=pd.merge(comp_Q[['gvkey','jdate','ROE']],ccm,how='left',on=['gvkey'])

# set link date bounds
ccm4=ccm3[(ccm3['jdate']>=ccm3['linkdt'])&(ccm3['jdate']<=ccm3['linkenddt'])].copy()
ccm4=ccm4[['gvkey','permno','jdate','ROE']]

# link comp and crsp using permno and jdate
ccm_ROE=pd.merge(crsp2[['permno', 'jdate']], ccm4, how='left', on=['permno', 'jdate'])
ccm_ROE=ccm_ROE.drop(['gvkey'],axis=1)
ccm_ROE['ROE'] = ccm_ROE['ROE'].fillna(method='ffill')

ccm_final = pd.merge(ccm_IA, ccm_ROE, how='left', on=['permno','jdate'])
ccm_final = ccm_final.loc[(ccm_final['jdate'] >= '1972-01-01') & (ccm_final['jdate'] <= '2018-11-01')]

ccm_final = ccm_final.loc[ccm_final['ROE'].notnull()]
ccm_final = ccm_final.loc[ccm_final['IA'].notnull()]
ccm_final = ccm_final.loc[ccm_final['lme'].notnull()]

#######################
# IA and ROE Return   #
#######################

# extract july data to calculate SIZE decile
ccm_july=ccm_final[ccm_final['ffmonth']==1].reset_index(drop=True)

# get size decile
# extract NYSE data to get breakpoints
ccm_NYSE = ccm_july.loc[ccm_july['exchcd'] == 1]
size_bins = ccm_NYSE.groupby('ffyear')['lme'].apply(lambda x: pd.qcut(x, 2, retbins=True)[1])
# Set interval boundary to -inf and inf
for i in size_bins.index:
    size_bins[i][0] = -np.inf
    size_bins[i][-1] = np.inf
# create size decile
ccm_july['size_decile'] = ccm_july.groupby('ffyear')['lme'].apply(lambda x: pd.cut(x, bins=size_bins[x.name], labels=np.arange(1, 3, 1)))

# same method to get IA decile
IA_bins = ccm_NYSE.groupby('ffyear')['IA'].apply(lambda x: pd.qcut(x, [0,0.3,0.7,1], retbins=True)[1])
# Set interval boundary to -inf and inf
for i in IA_bins.index:
    IA_bins[i][0] = -np.inf
    IA_bins[i][-1] = np.inf
ccm_july['IA_decile'] = ccm_july.groupby('ffyear')['IA'].apply(lambda x: pd.cut(x, bins=IA_bins[x.name], labels=np.arange(1, 4, 1)))

ccm_july=ccm_july[['permno','ffyear','size_decile','IA_decile']]

# merge decile back to ccm_final using FFyear. So same FFyear has the same decile.
ccm_decile = pd.merge(ccm_final,ccm_july,how='left',on=['permno','ffyear'])

ccm_decile['mdate'] = ccm_decile['jdate'].dt.year * 100 + ccm_decile['jdate'].dt.month
ccm_NYSE = ccm_decile.loc[ccm_decile['exchcd'] == 1]

# same method to get ROE decile
ROE_bins = ccm_NYSE.groupby('mdate')['ROE'].apply(lambda x: pd.qcut(x, [0,0.3,0.7,1], retbins=True)[1])

# Set interval boundary to -inf and inf
for i in ROE_bins.index:
    ROE_bins[i][0] = -np.inf
    ROE_bins[i][-1] = np.inf

ccm_decile['ROE_decile'] = ccm_decile.groupby('mdate')['ROE'].apply(lambda x: pd.cut(x, bins=ROE_bins[x.name], labels=np.arange(1, 4, 1)))


# calculate value weighted return for each decile and each month
# vw = lambda x: np.average(x, weights=ccm_decile.loc[x.index, "lme"])
# func_size = {'retadj':vw}
# func_IA = {'retadj':vw}
# func_ROE = {'retadj':vw}
#
# Ret_size = ccm_decile.groupby(['jdate', 'size_decile']).agg(func_size).reset_index()
# Ret_IA = ccm_decile.groupby(['jdate', 'IA_decile']).agg(func_IA).reset_index()
# Ret_ROE = ccm_decile.groupby(['jdate', 'ROE_decile']).agg(func_ROE).reset_index()

# Ret_size_pivot = Ret_size.pivot_table('retadj', ['jdate'], 'size_decile')
# Ret_size_pivot.columns = [1, 2]
# Ret_size_pivot=Ret_size_pivot.reset_index()
# Ret_size_pivot=Ret_size_pivot.sort_values(by=['jdate'])
#
# Ret_IA_pivot = Ret_IA.pivot_table('retadj', ['jdate'], 'IA_decile')
# Ret_IA_pivot.columns = [1, 2, 3]
# Ret_IA_pivot=Ret_IA_pivot.reset_index()
# Ret_IA_pivot=Ret_IA_pivot.sort_values(by=['jdate'])
#
# Ret_ROE_pivot = Ret_ROE.pivot_table('retadj', ['jdate'], 'ROE_decile')
# Ret_ROE_pivot.columns = [1, 2, 3]
# Ret_ROE_pivot=Ret_ROE_pivot.reset_index()
# Ret_ROE_pivot=Ret_ROE_pivot.sort_values(by=['jdate'])

ccm_decile = ccm_decile.loc[ccm_decile['size_decile'].notnull()]
ccm_decile = ccm_decile.loc[ccm_decile['IA_decile'].notnull()]
ccm_decile = ccm_decile.loc[ccm_decile['ROE_decile'].notnull()]
ccm_decile = ccm_decile.loc[ccm_decile['lme'].notnull()]

ccm_decile['Decile'] = ccm_decile['size_decile'].astype('str') + ccm_decile['IA_decile'].astype('str') + ccm_decile['ROE_decile'].astype('str')

vw = lambda x: np.average(x, weights=ccm_decile.loc[x.index, "lme"])
func_decile = {'retadj':vw}
Ret_Decile = ccm_decile.groupby(['jdate', 'Decile']).agg(func_decile).reset_index()
Ret_MKT = ccm_decile.groupby(['jdate']).agg(func_decile).reset_index()
Ret_MKT.columns = ['jdate', 'MKT']

Ret_Decile_pivot = Ret_Decile.pivot_table('retadj', ['jdate'], 'Decile')
Ret_Decile_pivot=Ret_Decile_pivot.reset_index()
Ret_Decile_pivot=Ret_Decile_pivot.sort_values(by=['jdate'])


Ret_Decile_pivot['ME'] = np.nanmean([Ret_Decile_pivot['111'], Ret_Decile_pivot['112'], Ret_Decile_pivot['113'], Ret_Decile_pivot['121'], Ret_Decile_pivot['122'], Ret_Decile_pivot['123'], Ret_Decile_pivot['131'], Ret_Decile_pivot['132'], Ret_Decile_pivot['133']], axis=0) - np.nanmean([Ret_Decile_pivot['211'], Ret_Decile_pivot['212'], Ret_Decile_pivot['213'], Ret_Decile_pivot['221'], Ret_Decile_pivot['222'], Ret_Decile_pivot['223'], Ret_Decile_pivot['231'], Ret_Decile_pivot['232'], Ret_Decile_pivot['233']], axis=0)
Ret_Decile_pivot['IA'] = np.nanmean([Ret_Decile_pivot['111'], Ret_Decile_pivot['112'], Ret_Decile_pivot['113'], Ret_Decile_pivot['211'], Ret_Decile_pivot['212'], Ret_Decile_pivot['213']], axis=0) - np.nanmean([Ret_Decile_pivot['131'], Ret_Decile_pivot['132'], Ret_Decile_pivot['133'], Ret_Decile_pivot['231'], Ret_Decile_pivot['232'], Ret_Decile_pivot['233']], axis=0)
Ret_Decile_pivot['ROE'] = np.nanmean([Ret_Decile_pivot['113'], Ret_Decile_pivot['123'], Ret_Decile_pivot['133'], Ret_Decile_pivot['213'], Ret_Decile_pivot['223'], Ret_Decile_pivot['233']], axis=0) - np.nanmean([Ret_Decile_pivot['111'], Ret_Decile_pivot['121'], Ret_Decile_pivot['131'], Ret_Decile_pivot['211'], Ret_Decile_pivot['221'], Ret_Decile_pivot['231']], axis=0)
Ret_Decile_pivot = Ret_Decile_pivot[['jdate', 'ME', 'IA', 'ROE']]

My_Qfactor = pd.merge(Ret_MKT, Ret_Decile_pivot, how='left', on=['jdate'])

#######################
#     Question 1      #
#######################

Q_factor = pd.read_excel('HXZ q-Factors (monthly 1967 to 2018).xlsx')
Q_factor = Q_factor.loc[66:621].reset_index(drop=True)
Q_factor['MKT'] = Q_factor['MKT-RF'] + Q_factor['RF']
Q_factor['MKT'].corr(My_Qfactor['MKT'])
Q_factor['ME'].corr(My_Qfactor['ME'])
Q_factor['IA'].corr(My_Qfactor['IA'])
Q_factor['ROE'].corr(My_Qfactor['ROE'])

FF = pd.read_csv('F-F_Research_Data_Factors.csv')
FF = FF.loc[555:1110].reset_index(drop=True)
FF.columns = ['date', 'Mkt-RF', 'SMB','HML','RF']
FF['date'] = pd.to_datetime(FF['date'],format='%Y%m')
FF['Mkt-RF'] = FF['Mkt-RF'].astype('float')
FF['RF'] = FF['RF'].astype('float')
FF['SMB'] = FF['SMB'].astype('float')
FF['HML'] = FF['HML'].astype('float')

FF_mom = pd.read_csv('F-F_Momentum_Factor_CSV.zip',skiprows=14,names=['date','mom'])[:1107].astype(float)
FF_mom['date'] = pd.to_datetime(FF_mom['date'],format='%Y%m')
FF_mom = FF_mom.loc[546:1101].reset_index(drop=True)

def stats_func(y,x):
    lm = sm.OLS(y, sm.add_constant(x))
    res = lm.fit()
    return res.params, res.tvalues, res.rsquared

# panel A
stats_func(Q_factor['ME'],sm.add_constant(Q_factor['MKT-RF']))
stats_func(Q_factor['ME'],sm.add_constant(FF[['Mkt-RF','SMB','HML']].astype('float')))
stats_func(Q_factor['ME'],sm.add_constant(pd.concat([FF[['Mkt-RF','SMB','HML']].astype('float'),FF_mom['mom']], axis=1)))
stats_func(Q_factor['IA'],sm.add_constant(Q_factor['MKT-RF']))
stats_func(Q_factor['IA'],sm.add_constant(FF[['Mkt-RF','SMB','HML']].astype('float')))
stats_func(Q_factor['IA'],sm.add_constant(pd.concat([FF[['Mkt-RF','SMB','HML']].astype('float'),FF_mom['mom']], axis=1)))
stats_func(Q_factor['ROE'],sm.add_constant(Q_factor['MKT-RF']))
stats_func(Q_factor['ROE'],sm.add_constant(FF[['Mkt-RF','SMB','HML']].astype('float')))
stats_func(Q_factor['ROE'],sm.add_constant(pd.concat([FF[['Mkt-RF','SMB','HML']].astype('float'),FF_mom['mom']], axis=1)))

# panel B
Q1_df = pd.concat([Q_factor[['ME','IA','ROE','MKT-RF']],FF[['SMB','HML']],FF_mom['mom']], axis=1)
Q1_panel_B = Q1_df.corr()



#######################
#     Question 2      #
#######################

My_Qfactor.columns = ['jdate', 'My_MKT', 'My_ME', 'My_IA', 'My_ROE']
My_Qfactor[['My_MKT', 'My_ME', 'My_IA', 'My_ROE']] *= 100
Q_factor.columns = ['Year', 'Month', 'Original_RF', 'Original_MKT-RF', 'Original_ME', 'Original_IA', 'Original_ROE', 'Original_MKT']
Q2_output = pd.concat([Q_factor, My_Qfactor], axis=1)
Q2_output['My_MKT-RF'] = Q2_output['My_MKT'] - Q2_output['Original_RF']
Q2_output = Q2_output[['jdate', 'Original_MKT-RF', 'My_MKT-RF', 'Original_ME', 'My_ME', 'Original_IA', 'My_IA', 'Original_ROE', 'My_ROE']]

Q2_panel_A = Q2_output[['Original_MKT-RF', 'My_MKT-RF', 'Original_ME', 'My_ME', 'Original_IA', 'My_IA', 'Original_ROE', 'My_ROE']].apply(lambda x:pd.Series({'Average':x.mean()*12, 'Volatility':x.std()*np.sqrt(12), 'Sharpe Ratio': np.sqrt(12)*x.mean()/x.std(), 'Skewness':x.skew(), 'Kurtosis':x.kurtosis()}))

Q2_output['Diff_MKT-RF'] = 100 * np.abs(Q2_output['My_MKT-RF'] - Q2_output['Original_MKT-RF'])
Q2_output['Diff_ME'] = 100 * np.abs(Q2_output['My_ME'] - Q2_output['Original_ME'])
Q2_output['Diff_IA'] = 100 * np.abs(Q2_output['My_IA'] - Q2_output['Original_IA'])
Q2_output['Diff_ROE'] = 100 * np.abs(Q2_output['My_ROE'] - Q2_output['Original_ROE'])

Q2_panel_B = Q2_output[['Diff_MKT-RF', 'Diff_ME', 'Diff_IA', 'Diff_ROE']].apply(lambda x:pd.Series({'Average':x.mean(), 'Volatility':x.std(), 'Min': x.min(), '25%': x.quantile(q=0.25), '50%': x.quantile(q=0.5), '75%': x.quantile(q=0.75),  'Max':x.max()}))


Q2_panel_C = pd.DataFrame(index=['correlation'],columns=['MKT-RF','ME', 'IA', 'ROE'])
Q2_panel_C['MKT-RF'] = Q2_output['My_MKT-RF'].corr(Q2_output['Original_MKT-RF'])
Q2_panel_C['ME'] = Q2_output['My_ME'].corr(Q2_output['Original_ME'])
Q2_panel_C['IA'] = Q2_output['My_IA'].corr(Q2_output['Original_IA'])
Q2_panel_C['ROE'] = Q2_output['My_ROE'].corr(Q2_output['Original_ROE'])


# Plot

plt.plot(Q2_output['Original_MKT-RF'])
plt.title('MKT-RF(Original)')
plt.show()
plt.plot(Q2_output['My_MKT-RF'])
plt.title('MKT-RF(Replicated)')
plt.show()
plt.plot(Q2_output['Diff_MKT-RF'])
plt.title('Difference between MKT-RF original and replication (abs diff)')
plt.show()

plt.plot(Q2_output['Original_ME'])
plt.title('ME(Original)')
plt.show()
plt.plot(Q2_output['My_ME'])
plt.title('ME(Replicated)')
plt.show()
plt.plot(Q2_output['Diff_ME'])
plt.title('Difference between ME original and replication (abs diff)')
plt.show()

plt.plot(Q2_output['Original_IA'])
plt.title('IA(Original)')
plt.show()
plt.plot(Q2_output['My_IA'])
plt.title('IA(Replicated)')
plt.show()
plt.plot(Q2_output['Diff_IA'])
plt.title('Difference between IA original and replication (abs diff)')
plt.show()

plt.plot(Q2_output['Original_ROE'])
plt.title('ROE(Original)')
plt.show()
plt.plot(Q2_output['My_ROE'])
plt.title('ROE(Replicated)')
plt.show()
plt.plot(Q2_output['Diff_ROE'])
plt.title('Difference between ROE original and replication (abs diff)')
plt.show()

#######################
#     Question 3      #
#######################

# panel A
stats_func(Q2_output['My_ME'],sm.add_constant(Q2_output['My_MKT-RF']))
stats_func(Q2_output['My_ME'],sm.add_constant(FF[['Mkt-RF','SMB','HML']].astype('float')))
stats_func(Q2_output['My_ME'],sm.add_constant(pd.concat([FF[['Mkt-RF','SMB','HML']].astype('float'),FF_mom['mom']], axis=1)))
stats_func(Q2_output['My_IA'],sm.add_constant(Q2_output['My_MKT-RF']))
stats_func(Q2_output['My_IA'],sm.add_constant(FF[['Mkt-RF','SMB','HML']].astype('float')))
stats_func(Q2_output['My_IA'],sm.add_constant(pd.concat([FF[['Mkt-RF','SMB','HML']].astype('float'),FF_mom['mom']], axis=1)))
stats_func(Q2_output['My_ROE'],sm.add_constant(Q2_output['My_MKT-RF']))
stats_func(Q2_output['My_ROE'],sm.add_constant(FF[['Mkt-RF','SMB','HML']].astype('float')))
stats_func(Q2_output['My_ROE'],sm.add_constant(pd.concat([FF[['Mkt-RF','SMB','HML']].astype('float'),FF_mom['mom']], axis=1)))

# panel B
Q1_df = pd.concat([Q_factor[['ME','IA','ROE','MKT-RF']],FF[['SMB','HML']],FF_mom['mom']], axis=1)
Q1_panel_B = Q1_df.corr()




