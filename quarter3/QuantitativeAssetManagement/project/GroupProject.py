import numpy as np
import pandas as pd
from pandas.tseries.offsets import *
import wrds
import matplotlib.pyplot as plt


###################
# Connect to WRDS #
###################
conn = wrds.Connection()

###################
# Compustat Data  #
###################
comp = conn.raw_sql("""
                    select gvkey, datadate, revt, cogs, at,
                    lt, pstkl, txditc, pstkrv, seq, pstk, 
                    ceq, mib, itcb, txdb, ni, dp, wcapch,
                    capx, sale, act, lct, che, dlc, txp,
                    dltt, mibt, pi, re, ebit, ib, csho
                    from comp.funda
                    where datafmt='STD'
                    and popsrc='D'
                    and consol='C'
                    and datadate >= '01/01/1950'
                    """)

comp['datadate']=pd.to_datetime(comp['datadate']) #convert datadate to date fmt
comp['year']=comp['datadate'].dt.year

# create Shareholders' equity
comp['she']=np.where(comp['seq'].isnull(), comp['ceq']+comp['pstk'], comp['seq'])
comp['she']=np.where(comp['she'].isnull(), comp['at']-comp['lt']-comp['mib'], comp['she'])
comp['she']=np.where(comp['she'].isnull(), comp['at']-comp['lt'], comp['she'])
comp['she']=np.where(comp['she'].isnull(), 0, comp['she'])

# Deferred taxes and investment tax credit
comp['dt']=np.where(comp['txditc'].isnull(), comp['itcb']+comp['txdb'], comp['txditc'])
comp['dt']=np.where(comp['dt'].isnull(), comp['itcb'], comp['dt'])
comp['dt']=np.where(comp['dt'].isnull(), comp['txdb'], comp['dt'])
comp['dt']=np.where(comp['dt'].isnull(), 0, comp['dt'])

# create preferrerd stock
comp['ps']=np.where(comp['pstkrv'].isnull(), comp['pstkl'], comp['pstkrv'])
comp['ps']=np.where(comp['ps'].isnull(),comp['pstk'], comp['ps'])
comp['ps']=np.where(comp['ps'].isnull(),0,comp['ps'])

# Add COMPUSTAT Pension
comp_pension = conn.raw_sql("""
                    select gvkey, datadate, prba
                    from comp.aco_pnfnda
                    where datafmt='STD'
                    and popsrc='D'
                    and consol='C'
                    and datadate >= '01/01/1950'
                    """)

comp_pension['datadate']=pd.to_datetime(comp['datadate']) #convert datadate to date fmt
comp = pd.merge(comp, comp_pension, how='left', on=['gvkey','datadate'])
comp['prba'] = comp['prba'].fillna(0)
# create book equity
comp['be']=comp['she']+comp['dt']-comp['ps']-comp['prba']
comp['be']=np.where(comp['be']>0, comp['be'], np.nan)
# same unit with market equity
comp['be']=comp['be']*1000
comp = comp.loc[comp['be'].notnull()]

comp=comp.sort_values(by=['gvkey','datadate'])

####### Profitability Variables #######
comp['gpoa']=(comp['revt']-comp['cogs'])/comp['at']
comp['roe']=comp['ib']/comp['be']
comp['roa']=comp['ib']/comp['at']
comp['wc']=comp['act']-comp['lct']-comp['che']+comp['dlc']+comp['txp']
comp['delta_wc']=comp['wc']-comp.groupby(['gvkey'])['wc'].shift(1)
comp['cfoa']=(comp['ni']+comp['dp']-comp['delta_wc']-comp['capx'])/comp['at']
comp['gmar']=(comp['revt']-comp['cogs'])/comp['sale']
comp['acc']=-(comp['delta_wc']-comp['dp'])/comp['at']

####### Growth Variables #######
comp['gp']=(comp['revt']-comp['cogs'])/comp['csho']
comp['delta_gpoa']=(comp['gp']-comp.groupby(['gvkey'])['gp'].shift(5))/comp.groupby(['gvkey'])['at'].shift(5)
comp['delta_roe']=(comp['ib']/comp['csho']-comp.groupby(['gvkey'])['ib'].shift(5)/comp.groupby(['gvkey'])['csho'].shift(5))/comp.groupby(['gvkey'])['be'].shift(5)
comp['delta_roa']=(comp['ib']/comp['csho']-comp.groupby(['gvkey'])['ib'].shift(5)/comp.groupby(['gvkey'])['csho'].shift(5))/comp.groupby(['gvkey'])['at'].shift(5)
comp['cf']=(comp['ib']+comp['dp']-comp['delta_wc']-comp['capx'])/comp['csho']
comp['delta_cfoa']=(comp['cf']-comp.groupby(['gvkey'])['cf'].shift(5))/comp.groupby(['gvkey'])['at'].shift(5)
comp['delta_gmar']=(comp['gp']/comp['csho']-comp.groupby(['gvkey'])['gp'].shift(5)/comp.groupby(['gvkey'])['csho'].shift(5))/comp.groupby(['gvkey'])['sale'].shift(5)

####### Safety Variables #######
comp['lev']=-(comp['dltt']+comp['dlc']+comp['mibt']+comp['pstk'])/comp['at']

comp=comp[['gvkey','datadate','year','be','at','lt','dlc','dltt','act','lct','ni','pi','wc','re','ebit','sale','roe','gpoa','roa','cfoa','gmar','acc','delta_gpoa','delta_roe','delta_roa','delta_cfoa','delta_gmar','lev','ib']]

###################
# CRSP Block      #
###################
crsp_m = conn.raw_sql("""
                      select a.permno, a.permco, a.date, b.shrcd, b.exchcd,
                      a.ret, a.shrout, a.prc
                      from crspa.msf as a
                      left join crspa.msenames as b
                      on a.permno=b.permno
                      and b.namedt<=a.date
                      and a.date<=b.nameendt
                      where a.date between '01/01/1950' and '12/31/2018'
                      and b.exchcd between 1 and 3
                      and b.shrcd between 10 and 11
                      """)

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
# sum of me across different permno belonging to same permco a given date
crsp_summe = crsp.groupby(['jdate','permco'])['me'].sum().reset_index()
# merge permco me to permno me
crsp=crsp.drop(['me'], axis=1)
crsp2=pd.merge(crsp, crsp_summe, how='left', on=['jdate','permco'])
# sort by permno and date
crsp2=crsp2.drop(['permco','date','ret'], axis=1)
crsp2=crsp2.sort_values(by=['permno','jdate'])

# t Jun to t+1 July is within a ffyear
crsp2['ffdate']=crsp2['jdate']+MonthEnd(-6)
crsp2['ffyear']=crsp2['ffdate'].dt.year
crsp2['ffmonth']=crsp2['ffdate'].dt.month

crsp2['year']=crsp2['jdate'].dt.year
crsp2['month']=crsp2['jdate'].dt.month

# lag market cap
crsp2['lme']=crsp2.groupby(['permno'])['me'].shift(1)

crsp_final=crsp2[['permno', 'jdate','year','month','ffyear','ffmonth','shrcd','exchcd','lme','me','retadj']]
crsp_final=crsp_final.sort_values(by=['permno','jdate'])

crsp_Dec=crsp_final[crsp_final['month']==12]

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

# to the year end before mering Compustat with CRSP (required by FF paper)
# "align accounting variables at the end of the firm’s fiscal year ending anywhere in calendar year t-1 to June of calendar year t"
ccm1=pd.merge(comp,ccm,how='left',on=['gvkey'])
ccm1['jdate']=ccm1['datadate']+YearEnd(0)

# set link date bounds
ccm2=ccm1[(ccm1['jdate']>=ccm1['linkdt'])&(ccm1['jdate']<=ccm1['linkenddt'])].copy()
# delete jdate, because jdate is not actual date any more
ccm2=ccm2.drop(['linktype','linkprim','linkdt','linkenddt','jdate'], axis=1)

# link comp and crsp Dec data using permno and year
ccm_final=pd.merge(crsp_Dec, ccm2, how='left', on=['permno', 'year'])
ccm_final=ccm_final.drop(['gvkey','datadate'], axis=1)

# Process missing Book equity data
#http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_historical_be_data.html
# be_missing = pd.read_fwf('DFF_BE_With_Nonindust.txt',header=None)
# be_missing = be_missing.drop([1,2],axis=1)
# be_missing.columns = ['permno'] + list(np.arange(1926, 2002))
# be_missing_melt = pd.melt(be_missing, id_vars=['permno'],value_vars=list(np.arange(1926, 2002)))
# be_missing_melt.columns = ['permno','ffyear','be']
# be_missing_melt = be_missing_melt.replace(-99.99,np.nan)
# be_missing_melt = be_missing_melt.dropna()
# ccm_final = pd.merge(ccm_final, be_missing_melt, how='left', on=['permno','ffyear'])
# ccm_final['be'] = np.where(ccm_final['be_x'].notnull(), ccm_final['be_x'], ccm_final['be_y'])

ccm_final=ccm_final.sort_values(by=['permno','jdate'])

####### Safety Variables #######
ccm_final['adjasset']=ccm_final['at']+0.1*(ccm_final['me']-ccm_final['be'])
ccm_final['tlta']=(ccm_final['dlc']+ccm_final['dltt'])/ccm_final['adjasset']
ccm_final['wcta']=(ccm_final['act']-ccm_final['lct'])/ccm_final['adjasset']
ccm_final['clca']=ccm_final['lct']-ccm_final['act']
ccm_final['oeneg']=(ccm_final['lt']>ccm_final['at']).astype(int)
ccm_final['nita']=ccm_final['ib']/ccm_final['at']
ccm_final['futl']=ccm_final['pi']/ccm_final['lt']
ccm_final['intwo']=(ccm_final['ib']<0).astype(int)
ccm_final['lib']=ccm_final.groupby(['permno'])['ib'].shift(1)
ccm_final['chin']=(ccm_final['ib']-ccm_final['lib'])/(np.abs(ccm_final['ib'])+np.abs(ccm_final['lib']))

# annual cpi, 2015 index equal to 100
CPI=pd.read_csv('CPALTT01USA661S.csv')
CPI['DATE']=pd.to_datetime(CPI['DATE'])
CPI['year']=CPI['DATE'].dt.year - 1
CPI.columns=['date','cpi','year']
CPI=CPI[['year','cpi']]
ccm_final=pd.merge(ccm_final,CPI,how='left',on=['year'])

ccm_final['O-Score']=-(-1.32-0.407*np.log(ccm_final['adjasset']/ccm_final['cpi'])+6.03*ccm_final['tlta']-1.43*ccm_final['wcta']+0.076*ccm_final['clca']-1.72*ccm_final['oeneg']-2.37*ccm_final['nita']-1.83*ccm_final['futl']+0.285*ccm_final['intwo']-0.521*ccm_final['chin'])
ccm_final['Z-Score']=(1.2*ccm_final['wc']+1.4*ccm_final['re']+3.3*ccm_final['ebit']+0.6*ccm_final['me']+ccm_final['sale'])/ccm_final['at']
ccm_final['evol']=ccm_final.groupby(['permno'])['roe'].rolling(5,min_periods=5).std().reset_index(0, drop=True)

ccm_final=ccm_final[['permno','jdate','gpoa','roe','roa','cfoa','gmar','acc','delta_gpoa','delta_roe','delta_roa','delta_cfoa','delta_gmar','lev','O-Score','Z-Score','evol']]

####### Quality Z scores #######
ccm_final['r_gpoa']=ccm_final.groupby(['jdate'])['gpoa'].rank()
ccm_final['r_roe']=ccm_final.groupby(['jdate'])['roe'].rank()
ccm_final['r_roa']=ccm_final.groupby(['jdate'])['roa'].rank()
ccm_final['r_cfoa']=ccm_final.groupby(['jdate'])['cfoa'].rank()
ccm_final['r_gmar']=ccm_final.groupby(['jdate'])['gmar'].rank()
ccm_final['r_acc']=ccm_final.groupby(['jdate'])['acc'].rank()
ccm_final['r_dgpoa']=ccm_final.groupby(['jdate'])['delta_gpoa'].rank()
ccm_final['r_droe']=ccm_final.groupby(['jdate'])['delta_roe'].rank()
ccm_final['r_droa']=ccm_final.groupby(['jdate'])['delta_roa'].rank()
ccm_final['r_dcfoa']=ccm_final.groupby(['jdate'])['delta_cfoa'].rank()
ccm_final['r_dgmar']=ccm_final.groupby(['jdate'])['delta_gmar'].rank()
ccm_final['r_lev']=ccm_final.groupby(['jdate'])['lev'].rank()
ccm_final['r_o']=ccm_final.groupby(['jdate'])['O-Score'].rank()
ccm_final['r_z']=ccm_final.groupby(['jdate'])['Z-Score'].rank()
ccm_final['r_evol']=ccm_final.groupby(['jdate'])['evol'].rank()


ccm_final['z_gpoa']=(ccm_final['r_gpoa']-ccm_final.groupby(['jdate'])['r_gpoa'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_gpoa'].transform(np.std)
ccm_final['z_roe']=(ccm_final['r_roe']-ccm_final.groupby(['jdate'])['r_roe'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_roe'].transform(np.std)
ccm_final['z_roa']=(ccm_final['r_roa']-ccm_final.groupby(['jdate'])['r_roa'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_roa'].transform(np.std)
ccm_final['z_cfoa']=(ccm_final['r_cfoa']-ccm_final.groupby(['jdate'])['r_cfoa'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_cfoa'].transform(np.std)
ccm_final['z_gmar']=(ccm_final['r_gmar']-ccm_final.groupby(['jdate'])['r_gmar'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_gmar'].transform(np.std)
ccm_final['z_acc']=(ccm_final['r_acc']-ccm_final.groupby(['jdate'])['r_acc'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_acc'].transform(np.std)
ccm_final['z_dgpoa']=(ccm_final['r_dgpoa']-ccm_final.groupby(['jdate'])['r_dgpoa'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_dgpoa'].transform(np.std)
ccm_final['z_droe']=(ccm_final['r_droe']-ccm_final.groupby(['jdate'])['r_droe'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_droe'].transform(np.std)
ccm_final['z_droa']=(ccm_final['r_droa']-ccm_final.groupby(['jdate'])['r_droa'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_droa'].transform(np.std)
ccm_final['z_dcfoa']=(ccm_final['r_dcfoa']-ccm_final.groupby(['jdate'])['r_dcfoa'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_dcfoa'].transform(np.std)
ccm_final['z_dgmar']=(ccm_final['r_dgmar']-ccm_final.groupby(['jdate'])['r_dgmar'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_dgmar'].transform(np.std)
ccm_final['z_lev']=(ccm_final['r_lev']-ccm_final.groupby(['jdate'])['r_lev'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_lev'].transform(np.std)
ccm_final['z_o']=(ccm_final['r_o']-ccm_final.groupby(['jdate'])['r_o'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_o'].transform(np.std)
ccm_final['z_z']=(ccm_final['r_z']-ccm_final.groupby(['jdate'])['r_z'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_z'].transform(np.std)
ccm_final['z_evol']=(ccm_final['r_evol']-ccm_final.groupby(['jdate'])['r_evol'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_evol'].transform(np.std)

ccm_final['profitability']=ccm_final['z_gpoa']+ccm_final['z_roe']+ccm_final['z_roa']+ccm_final['z_cfoa']+ccm_final['z_gmar']+ccm_final['z_acc']
ccm_final['growth']=ccm_final['z_dgpoa']+ccm_final['z_droe']+ccm_final['z_droa']+ccm_final['z_dcfoa']+ccm_final['z_dgmar']
# negative because low volatility should have higher z scores
ccm_final['safety']=-(ccm_final['z_lev']+ccm_final['z_o']+ccm_final['z_z']+ccm_final['z_evol'])

ccm_final['r_profit']=ccm_final.groupby(['jdate'])['profitability'].rank()
ccm_final['r_growth']=ccm_final.groupby(['jdate'])['growth'].rank()
ccm_final['r_safety']=ccm_final.groupby(['jdate'])['safety'].rank()

ccm_final['z_profit']=(ccm_final['r_profit']-ccm_final.groupby(['jdate'])['r_profit'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_profit'].transform(np.std)
ccm_final['z_growth']=(ccm_final['r_growth']-ccm_final.groupby(['jdate'])['r_growth'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_growth'].transform(np.std)
ccm_final['z_safety']=(ccm_final['r_safety']-ccm_final.groupby(['jdate'])['r_safety'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_safety'].transform(np.std)

ccm_final['quality']=ccm_final['z_profit']+ccm_final['z_growth']+ccm_final['z_safety']

ccm_final['r_quality']=ccm_final.groupby(['jdate'])['quality'].rank()
ccm_final['z_quality']=(ccm_final['r_quality']-ccm_final.groupby(['jdate'])['r_quality'].transform(np.mean))/ccm_final.groupby(['jdate'])['r_quality'].transform(np.std)

ccm_quality=ccm_final[['permno','jdate','z_quality']].copy()

# Merge yearly quality with crsp data using fiscal year.
# add 1 because we have year t-1 data on July of year t
ccm_quality['ffyear']=ccm_quality['jdate'].dt.year + 1
ccm_quality=ccm_quality.drop(['jdate'], axis=1)
ccm_all=pd.merge(crsp_final,ccm_quality,how='left',on=['permno','ffyear'])

# extract july data to calculate decile
ccm_july=ccm_all[ccm_all['month']==7].reset_index(drop=True)
ccm_july=ccm_july[ccm_july['z_quality'].notnull()]
# extract NYSE data to get breakpoints
ccm_NYSE = ccm_july.loc[ccm_july['exchcd'] == 1]
quality_bins = ccm_NYSE.groupby('year')['z_quality'].apply(lambda x: pd.qcut(x, 10, retbins=True)[1])
# Set interval boundary to -inf and inf
for i in quality_bins.index:
    quality_bins[i][0] = -np.inf
    quality_bins[i][-1] = np.inf
# create quality decile
ccm_july['quality_decile'] = ccm_july.groupby('year')['z_quality'].apply(lambda x: pd.cut(x, bins=quality_bins[x.name], labels=np.arange(1, 11, 1)))
ccm_july=ccm_july[['permno','ffyear','quality_decile']]

# merge decile back to ccm_all using FFyear. So same FFyear has the same decile.
ccm_decile = pd.merge(ccm_all,ccm_july,how='left',on=['permno','ffyear'])

# calculate value weighted return for each decile and each month
vw = lambda x: np.average(x, weights=ccm_decile.loc[x.index, "lme"])
func_quality = {'retadj':vw}
Ret_quality = ccm_decile.groupby(['jdate', 'quality_decile']).agg(func_quality).reset_index()
# get pivot table
Ret_quality_pivot = Ret_quality.pivot_table('retadj', ['jdate'], 'quality_decile')
Ret_quality_pivot.columns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
Ret_quality_pivot=Ret_quality_pivot.reset_index()
Ret_quality_pivot=Ret_quality_pivot.sort_values(by=['jdate'])
Ret_quality_pivot[11]=Ret_quality_pivot[10]-Ret_quality_pivot[1]


###################
# Official Data   #
###################

AQR_10port = pd.read_excel('Quality Minus Junk 10 QualitySorted Portfolios Monthly.xlsx', skiprows=18)
AQR_10port = AQR_10port[AQR_10port.columns[0:11]]
AQR_10port.columns = ['date', 1, 2, 3, 4, 5, 6, 7 ,8 ,9, 10]
AQR_10port[11] = AQR_10port[10] - AQR_10port[1]
AQR_10port = AQR_10port.loc[0:737]

AQR_QMJ = pd.read_excel('Quality Minus Junk Factors Monthly.xlsx', skiprows=18)
AQR_QMJ = AQR_QMJ[['DATE', 'USA']]
AQR_QMJ = AQR_QMJ.loc[0:737]


FF = pd.read_csv('/Users/leonard/Desktop/MFE-QAM/PS/F-F_Research_Data_Factors.csv')
FF = FF.loc[507:1112].reset_index(drop=True)
FF.columns = ['date', 'Mkt-RF', 'SMB','HML','RF']
FF[['Mkt-RF','SMB','HML','RF']] = FF[['Mkt-RF','SMB','HML','RF']].astype('float')
FF['Mkt'] = FF['Mkt-RF'] + FF['RF']
FF['date'] = pd.to_datetime(FF['date'],format='%Y%m')
FF['date'] = FF['date'] + MonthEnd(0)


### Cumulative Return
Quality_cumret = np.exp(np.cumsum(np.log(Ret_quality_pivot[11]+1)))
Mkt_cumret = np.exp(np.cumsum(np.log(FF['Mkt-RF']/100+1)))
Qua = np.exp(np.cumsum(np.log(AQR_10port.loc[132:737,11].reset_index(drop=True)+1)))

plt.plot(Qua,'r',label='Quality_AQR')
plt.plot(Quality_cumret,'k',label='Quality')
plt.plot(Mkt_cumret,'b',label='Market')
plt.legend()
plt.show()

### Correlations
Ret_quality_pivot[10].corr(AQR_10port.loc[132:737,10])

np.mean(Ret_quality_pivot[11])
np.mean(AQR_10port.loc[132:737,11])
