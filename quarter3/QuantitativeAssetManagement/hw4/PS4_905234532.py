import pandas as pd
import numpy as np
import wrds
from pandas.tseries.offsets import *
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import skew, kurtosis


conn = wrds.Connection()

comp = conn.raw_sql("""
                    select gvkey, datadate, at, lt, pstkl, txditc,
                    pstkrv, seq, pstk, ceq, mib, itcb, txdb
                    from comp.funda
                    where indfmt='INDL'
                    and datafmt='STD'
                    and popsrc='D'
                    and consol='C'
                    and datadate >= '01/01/1970'
                    """)

comp['datadate']=pd.to_datetime(comp['datadate'])
comp['year']=comp['datadate'].dt.year

comp['she']=np.where(comp['seq'].isnull(), comp['ceq']+comp['pstk'], comp['seq'])
comp['she']=np.where(comp['she'].isnull(), comp['at']-comp['lt']-comp['mib'], comp['she'])
comp['she']=np.where(comp['she'].isnull(), comp['at']-comp['lt'], comp['she'])
comp['she']=np.where(comp['she'].isnull(), 0, comp['she'])


comp['dt']=np.where(comp['txditc'].isnull(), comp['itcb']+comp['txdb'], comp['txditc'])
comp['dt']=np.where(comp['dt'].isnull(), comp['itcb'], comp['dt'])
comp['dt']=np.where(comp['dt'].isnull(), comp['txdb'], comp['dt'])
comp['dt']=np.where(comp['dt'].isnull(), 0, comp['dt'])


comp['ps']=np.where(comp['pstkrv'].isnull(), comp['pstkl'], comp['pstkrv'])
comp['ps']=np.where(comp['ps'].isnull(),comp['pstk'], comp['ps'])
comp['ps']=np.where(comp['ps'].isnull(),0,comp['ps'])


comp_pension = conn.raw_sql("""
                    select gvkey, datadate, prba
                    from comp.aco_pnfnda
                    where indfmt='INDL' 
                    and datafmt='STD'
                    and popsrc='D'
                    and consol='C'
                    and datadate >= '01/01/1959'
                    """)

comp_pension['datadate']=pd.to_datetime(comp['datadate']) #convert datadate to date fmt
comp = pd.merge(comp, comp_pension, how='left', on=['gvkey','datadate'])
comp['prba'] = comp['prba'].fillna(0)
# create book equity
comp['be']=comp['she']+comp['dt']-comp['ps']-comp['prba']
comp['be']=np.where(comp['be']>0, comp['be'], np.nan)
comp = comp.loc[comp['be'].notnull()]

comp=comp[['gvkey','datadate','year','be']]


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

# change variable format to int
crsp_m[['permco','permno','shrcd','exchcd']]=crsp_m[['permco','permno','shrcd','exchcd']].astype(int)

# Line up date to be end of month
crsp_m['date']=pd.to_datetime(crsp_m['date'])
crsp_m['jdate']=crsp_m['date']+MonthEnd(0)

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


crsp_summe = crsp.groupby(['jdate','permco'])['me'].sum().reset_index()
# merge permco me to permno me
crsp=crsp.drop(['me'], axis=1)
crsp2=pd.merge(crsp, crsp_summe, how='left', on=['jdate','permco'])
# sort by permno and date
crsp2=crsp2.drop(['permco','date','ret'], axis=1)
crsp2=crsp2.sort_values(by=['permno','jdate'])

# keep December market cap
crsp2['year']=crsp2['jdate'].dt.year
crsp2['month']=crsp2['jdate'].dt.month
decme=crsp2[crsp2['month']==12]
decme=decme[['permno','jdate','me','year']].rename(columns={'me':'dec_me'})

# t Jun to t+1 July is within a ffyear
crsp2['ffdate']=crsp2['jdate']+MonthEnd(-6)
crsp2['ffyear']=crsp2['ffdate'].dt.year
crsp2['ffmonth']=crsp2['ffdate'].dt.month

# lag market cap
crsp2['lme']=crsp2.groupby(['permno'])['me'].shift(1)
crsp2=crsp2.drop(['me'], axis=1)

# decme shift(1)
decme['ffyear']=decme['year']+1
decme=decme[['permno','ffyear','dec_me']]


crsp_final = pd.merge(crsp2, decme, how='left', on=['permno','ffyear'])
crsp_final=crsp_final[['permno', 'jdate','ffyear','ffmonth','shrcd','exchcd','lme','dec_me','retadj']]
crsp_final=crsp_final.sort_values(by=['permno','jdate'])


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


ccm1=pd.merge(comp[['gvkey','datadate','be']],ccm,how='left',on=['gvkey'])
ccm1['yearend']=ccm1['datadate']+YearEnd(0)
ccm1['jdate']=ccm1['yearend']+MonthEnd(6)

# set link date bounds
ccm2=ccm1[(ccm1['jdate']>=ccm1['linkdt'])&(ccm1['jdate']<=ccm1['linkenddt'])].copy()
# create ffyear and delete jdate, because jdate is not actual date any more
ccm2['ffyear']=ccm2['jdate'].dt.year
ccm2=ccm2[['gvkey','permno','ffyear','be']]

# link comp and crsp using permno and ffyear
ccm_final=pd.merge(crsp_final, ccm2, how='left', on=['permno', 'ffyear'])
ccm_final=ccm_final.drop(['gvkey'],axis=1)
ccm_final['year']=ccm_final['jdate'].dt.year
ccm_final['month']=ccm_final['jdate'].dt.month




# extract july data to calculate decile
ccm_july=ccm_final[ccm_final['month']==7].reset_index(drop=True)
ccm_july['bm']=ccm_july['be']*1000/ccm_july['dec_me']
ccm_july=ccm_july[ccm_july['bm'].notnull()]


ccm_NYSE = ccm_july.loc[ccm_july['exchcd'] == 1]
size_bins = ccm_NYSE.groupby('year')['lme'].apply(lambda x: pd.qcut(x, 10, retbins=True)[1])

ccm_july['size_decile'] = ccm_july[['lme','year']].apply(lambda x: np.digitize(x['lme'],size_bins[x['year']]),axis=1).replace({11:10,0:1})

bm_bins = ccm_NYSE.groupby('year')['bm'].apply(lambda x: pd.qcut(x, 10, retbins=True)[1])

ccm_july['bm_decile'] = ccm_july[['bm','year']].apply(lambda x: np.digitize(x['bm'],bm_bins[x['year']]),axis=1).replace({11:10,0:1})
ccm_july=ccm_july[['permno','ffyear','size_decile','bm_decile']]

ccm_decile = pd.merge(ccm_final,ccm_july,how='left',on=['permno','ffyear'])

# calculate value weighted return for each decile and each month
vw = lambda x: np.average(x, weights=ccm_decile.loc[x.index, "lme"])
func_size = {'retadj':vw}
func_bm = {'retadj':vw}
Ret_size = ccm_decile.groupby(['jdate', 'size_decile']).agg(func_size).reset_index()
# get pivot table
Ret_size_pivot = Ret_size.pivot_table('retadj', ['jdate'], 'size_decile')
Ret_size_pivot.columns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
Ret_size_pivot=Ret_size_pivot.reset_index()
Ret_size_pivot=Ret_size_pivot.sort_values(by=['jdate'])
# Jan 1973 to Dec 2018
Ret_size_pivot=Ret_size_pivot[18:570].reset_index(drop=True)
# decimal to percentage
Ret_size_pivot[[1,2,3,4,5,6,7,8,9,10]] = Ret_size_pivot[[1,2,3,4,5,6,7,8,9,10]]*100
# 11 is Long-Short portfolio
Ret_size_pivot[11] = Ret_size_pivot[1] - Ret_size_pivot[10]


Ret_bm = ccm_decile.groupby(['jdate', 'bm_decile']).agg(func_bm).reset_index()
Ret_bm_pivot = Ret_bm.pivot_table('retadj', ['jdate'], 'bm_decile')
Ret_bm_pivot.columns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
Ret_bm_pivot=Ret_bm_pivot.reset_index()
Ret_bm_pivot=Ret_bm_pivot.sort_values(by=['jdate'])
# Jan 1973 to Dec 2018
Ret_bm_pivot=Ret_bm_pivot[18:570].reset_index(drop=True)
Ret_bm_pivot[[1,2,3,4,5,6,7,8,9,10]] = Ret_bm_pivot[[1,2,3,4,5,6,7,8,9,10]]*100
Ret_bm_pivot[11] = Ret_bm_pivot[10] - Ret_bm_pivot[1]



ccm_july_6factor=ccm_final[ccm_final['month']==7].reset_index(drop=True)
ccm_july_6factor['bm']=ccm_july_6factor['be']*1000/ccm_july_6factor['dec_me']
ccm_july_6factor=ccm_july_6factor[ccm_july_6factor['bm'].notnull()]

# Divide size to 2 parts
size_bins = ccm_NYSE.groupby('year')['lme'].apply(lambda x: pd.qcut(x, 2, retbins=True)[1])
# Set interval boundary to -inf and inf
for i in size_bins.index:
    size_bins[i][0] = -np.inf
    size_bins[i][-1] = np.inf

ccm_july_6factor['size_decile'] = ccm_july_6factor.groupby('year')['lme'].apply(lambda x: pd.cut(x, bins=size_bins[x.name], labels=np.arange(1, 3, 1)))

bm_bins = ccm_NYSE.groupby('year')['bm'].apply(lambda x: pd.qcut(x, [0,0.3,0.7,1], retbins=True)[1])
# Set interval boundary to -inf and inf
for i in bm_bins.index:
    bm_bins[i][0] = -np.inf
    bm_bins[i][-1] = np.inf
ccm_july_6factor['bm_decile'] = ccm_july_6factor.groupby('year')['bm'].apply(lambda x: pd.cut(x, bins=bm_bins[x.name], labels=np.arange(1, 4, 1)))
# Replace decile number with string
map_dict_size = {1:'Small', 2:'Big'}
map_dict_bm = {1:'Growth', 2:'Neutral', 3:'Value'}
ccm_july_6factor['size_decile']=ccm_july_6factor['size_decile'].map(map_dict_size)
ccm_july_6factor['bm_decile']=ccm_july_6factor['bm_decile'].map(map_dict_bm)
ccm_july_6factor['decile']=ccm_july_6factor['size_decile']+' '+ccm_july_6factor['bm_decile']
ccm_july_6factor=ccm_july_6factor[['permno','ffyear','decile']]

# merge decile back to all data using FFyear
ccm_6factor = pd.merge(ccm_final,ccm_july_6factor,how='left',on=['permno','ffyear'])
# calculate value weighted return for each decile and each month
func_6factor = {'retadj':vw}
Ret_6factor = ccm_6factor.groupby(['jdate', 'decile']).agg(func_6factor).reset_index()
Ret_6factor_pivot = Ret_6factor.pivot_table('retadj', ['jdate'], 'decile')
Ret_6factor_pivot=Ret_6factor_pivot.reset_index()
Ret_6factor_pivot=Ret_6factor_pivot.sort_values(by=['jdate'])
# Jan 1973 to Dec 2018
Ret_6factor_pivot=Ret_6factor_pivot[18:570].reset_index(drop=True)
Ret_6factor_pivot[['Big Growth','Big Neutral','Big Value','Small Growth','Small Neutral','Small Value']] = Ret_6factor_pivot[['Big Growth','Big Neutral','Big Value','Small Growth','Small Neutral','Small Value']]*100
Ret_6factor_pivot=Ret_6factor_pivot[['jdate','Small Growth','Small Neutral','Small Value','Big Growth','Big Neutral','Big Value']]
Ret_6factor_pivot['SMB'] = 1/3*(Ret_6factor_pivot['Small Value']+Ret_6factor_pivot['Small Neutral']+Ret_6factor_pivot['Small Growth']) - 1/3*(Ret_6factor_pivot['Big Value']+Ret_6factor_pivot['Big Neutral']+Ret_6factor_pivot['Big Growth'])
Ret_6factor_pivot['HML'] = 1/2*(Ret_6factor_pivot['Small Value']+Ret_6factor_pivot['Big Value']) - 1/2*(Ret_6factor_pivot['Small Growth']+Ret_6factor_pivot['Big Growth'])


Q1_output = pd.merge(Ret_bm, Ret_size, how='inner', left_on=['jdate','bm_decile'], right_on=['jdate','size_decile'])
Q1_output = Q1_output.drop(['size_decile'],axis=1)
Q1_output = pd.merge(Q1_output,Ret_6factor_pivot[['jdate','SMB','HML']],how='left',on=['jdate'])
Q1_output.columns = ['jdate', 'port', 'BtM_Ret','Size_Ret','SMB_Ret','HML_Ret']
Q1_output['Year'] = Q1_output['jdate'].dt.year
Q1_output['Month'] = Q1_output['jdate'].dt.month
Q1_output = Q1_output[Q1_output['Year']>=1973]
Q1_output = Q1_output.drop(['jdate'],axis=1).reset_index(drop=True)
Q1_output = Q1_output[['Year','Month','port','Size_Ret','BtM_Ret','HML_Ret','SMB_Ret']]


FF = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/hw4/F-F_Research_Data_Factors.CSV',skiprows=3)
FF = FF.loc[558:1109].reset_index(drop=True)
FF.columns = ['date', 'Mkt-RF', 'SMB','HML','RF']
FF[['Mkt-RF','SMB','HML','RF']] = FF[['Mkt-RF','SMB','HML','RF']].astype('float')
FF['date'] = pd.to_datetime(FF['date'],format='%Y%m')
FF['date'] = FF['date'] + MonthEnd(0)

FF_size = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/hw4/Portfolios_Formed_on_ME_CSV.zip',skiprows=12)
FF_size_decile=FF_size[['Unnamed: 0','Lo 10','Dec 2','Dec 3','Dec 4','Dec 5','Dec 6','Dec 7','Dec 8','Dec 9','Hi 10']]
FF_size_decile.columns=['date',1,2,3,4,5,6,7,8,9,10]
# Jan 1973 to Dec 2018
FF_size_decile=FF_size_decile[558:1110].reset_index(drop=True)
FF_size_decile[[1,2,3,4,5,6,7,8,9,10]]=FF_size_decile[[1,2,3,4,5,6,7,8,9,10]].astype('float')
FF_size_decile[11] = FF_size_decile[1] - FF_size_decile[10]



corr_size = pd.DataFrame(FF_size_decile.iloc[:,1:12].corrwith(Ret_size_pivot.iloc[:,1:12]),columns=['correlation']).transpose()

Ret_size_pivot_ex = pd.merge(Ret_size_pivot,FF[['date','RF']],how='left',left_on=['jdate'],right_on=['date'])
Ret_size_pivot_ex[[1,2,3,4,5,6,7,8,9,10]] = Ret_size_pivot_ex[[1,2,3,4,5,6,7,8,9,10]].apply(lambda x: x-Ret_size_pivot_ex['RF'])
Q2_output = Ret_size_pivot_ex[[1,2,3,4,5,6,7,8,9,10,11]].apply(lambda x:pd.Series({'Annualized Mean':x.mean()*12, 'Annualized Volatility':x.std()*np.sqrt(12), 'Sharpe Ratio': np.sqrt(12)*x.mean()/x.std(), 'Skewness':x.skew()}))
Q2_output = Q2_output.append(corr_size)
Q2_output.columns = [1,2,3,4,5,6,7,8,9,10,'Long-Short']


FF_bm = pd.read_csv('/Users/huanyu/Desktop/QuantitativeAssetManagement/hw4/Portfolios_Formed_on_BE-ME_CSV.zip',skiprows=23)
FF_bm_decile=FF_bm[['Unnamed: 0','Lo 10','Dec 2','Dec 3','Dec 4','Dec 5','Dec 6','Dec 7','Dec 8','Dec 9','Hi 10']]
FF_bm_decile.columns=['date',1,2,3,4,5,6,7,8,9,10]
# Jan 1973 to Dec 2018
FF_bm_decile=FF_bm_decile[558:1110].reset_index(drop=True)
FF_bm_decile[[1,2,3,4,5,6,7,8,9,10]]=FF_bm_decile[[1,2,3,4,5,6,7,8,9,10]].astype('float')
FF_bm_decile[11] = FF_bm_decile[10] - FF_bm_decile[1]



corr_bm = pd.DataFrame(FF_bm_decile.iloc[:,1:12].corrwith(Ret_bm_pivot.iloc[:,1:12]),columns=['correlation']).transpose()

Ret_bm_pivot_ex = pd.merge(Ret_bm_pivot,FF[['date','RF']],how='left',left_on=['jdate'],right_on=['date'])
Ret_bm_pivot_ex[[1,2,3,4,5,6,7,8,9,10]] = Ret_bm_pivot_ex[[1,2,3,4,5,6,7,8,9,10]].apply(lambda x: x-Ret_bm_pivot_ex['RF'])
Q3_output = Ret_bm_pivot_ex[[1,2,3,4,5,6,7,8,9,10,11]].apply(lambda x:pd.Series({'Annualized Mean':x.mean()*12, 'Annualized Volatility':x.std()*np.sqrt(12), 'Sharpe Ratio': np.sqrt(12)*x.mean()/x.std(), 'Skewness':x.skew()}))
Q3_output = Q3_output.append(corr_bm)
Q3_output.columns = [1,2,3,4,5,6,7,8,9,10,'Long-Short']

Size_return = Ret_size_pivot[['jdate',11]][444:].reset_index(drop=True)
Value_return = Ret_bm_pivot[['jdate',11]][444:].reset_index(drop=True)
Size_return.mean()
Value_return.mean()
Size_cumret = np.exp(np.cumsum(np.log(Size_return[11]/100+1)))
Value_cumret = np.exp(np.cumsum(np.log(Value_return[11]/100+1)))

plt.plot(Size_cumret,'r',label='Size')
plt.plot(Value_cumret,'b',label='Value')
plt.legend()
plt.show()
plt.close()

corr_SMBHML = pd.DataFrame(FF[['SMB','HML']].corrwith(Ret_6factor_pivot[['SMB','HML']]),columns=['correlation']).transpose()


Q5_output = Ret_6factor_pivot[['SMB','HML']].apply(lambda x:pd.Series({'Annualized Mean':x.mean()*12, 'Annualized Volatility':x.std()*np.sqrt(12), 'Sharpe Ratio': np.sqrt(12)*x.mean()/x.std(), 'Skewness':x.skew()}))
Q5_output = Q5_output.append(corr_SMBHML)

#######################




