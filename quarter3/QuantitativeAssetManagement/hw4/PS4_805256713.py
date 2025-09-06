#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 14:14:23 2019

@author: Wesley
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

Crsp = pd.read_csv('~/Desktop/QAM/PS4/crsp.csv')
Compustat = pd.read_csv('~/Desktop/QAM/PS4/compustat.csv')

Compustat['Year'] = Compustat['datadate']//10000
Compustat['Month'] = Compustat['datadate']//100%100

Compustat = Compustat.loc[(Compustat.datafmt=="STD") & (Compustat.indfmt=="INDL")]

Pension = pd.read_csv('~/Desktop/QAM/PS4/pension.csv')
Pension['Year'] = Pension['datadate']//10000

# merge pension and compustat
Compustat = pd.merge(Compustat, Pension, on=['gvkey', 'Year'], how='outer')

# Shareholder's equity
Compustat['she'] = np.where(Compustat['seq'].isnull(), Compustat['ceq']+Compustat['pstk'], Compustat['seq'])
Compustat['she'] = np.where(Compustat['she'].isnull(), Compustat['at']-Compustat['lt']-Compustat['mib'], Compustat['she'])
Compustat['she'] = np.where(Compustat['she'].isnull(), Compustat['at']-Compustat['lt'], Compustat['she'])

# Deferred taxes and investment tax credit
Compustat['dt'] = np.where(Compustat['txditc'].isnull(), Compustat['itcb']+Compustat['txdb'], Compustat['txditc'])
Compustat['dt'] = np.where(Compustat['dt'].isnull(), Compustat['itcb'], Compustat['dt'])
Compustat['dt'] = np.where(Compustat['dt'].isnull(), Compustat['txdb'], Compustat['dt'])

# preferred stock
Compustat['ps']=np.where(Compustat['pstkrv'].isnull(), Compustat['pstkl'], Compustat['pstkrv'])
Compustat['ps']=np.where(Compustat['ps'].isnull(), Compustat['pstk'], Compustat['ps'])

# Book equity
Compustat['be'] = Compustat['she']-Compustat['ps']+Compustat['dt']-Compustat['prba']
Compustat['be'] = np.where(Compustat['be']>0, Compustat['be'], np.nan)

# number of years in Compustat
Compustat = Compustat.sort_values(by=['gvkey','Year','Month'])
Compustat = Compustat[['gvkey','Year','Month','be','datadate_x']]

# Following for dealing with CRSP
Crsp = Crsp.loc[~Crsp['RET'].isin(['C', 'B'])]
Crsp = Crsp.loc[~Crsp['DLRET'].isin(['A', 'S', 'T', 'P'])]
Crsp = Crsp.loc[Crsp['SHRCD'].isin([10, 11])]
Crsp = Crsp.loc[Crsp['EXCHCD'].isin([1, 2, 3])]
Crsp['RET'] = Crsp['RET'].fillna(0)
Crsp['DLRET'] = Crsp['DLRET'].fillna(0)
Crsp['RET'] = Crsp['RET'].astype('float')
Crsp['DLRET'] = Crsp['DLRET'].astype('float')
# cum-dividend returns
Crsp["RET"] = (1+Crsp.RET) * (1+Crsp.DLRET) - 1
Crsp['Mkt_Cap'] = Crsp['PRC'].abs() * Crsp['SHROUT']/1000
Crsp['Year'] = Crsp['date']//10000
Crsp['Month'] = Crsp['date']//100%100
# aggregate 
Crsp_monthly = Crsp.groupby(['Year','Month','PERMCO'])['Mkt_Cap'].sum()



# linked table linking and cleaning it
Link = pd.read_csv('~/Desktop/QAM/PS4/link.csv')
Link = Link.loc[(Link.LINKPRIM=='P') | (Link.LINKPRIM=='C')]
Link = Link[['gvkey','LPERMCO','LPERMNO','LINKDT','LINKENDDT']]
Link.LINKDT = pd.to_datetime(Link.LINKDT, format='%Y%m%d')
Link.LINKENDDT = Link.LINKENDDT.replace('E', "20181231")
Link.LINKENDDT = pd.to_datetime(Link.LINKENDDT, format='%Y%m%d')
# merge Compustate and Link
combo = pd.merge(Compustat, Link, on='gvkey', how='outer')
combo.datadate_x = pd.to_datetime(combo.datadate_x, format='%Y%m%d')
combo = combo.loc[(combo.datadate_x >= combo.LINKDT) & (combo.datadate_x <= combo.LINKENDDT)]

# get the non financial firm PERMCOS
Non_financial = combo.LPERMCO.unique()

final_merge = pd.merge(Crsp, combo, left_on=['PERMCO','Year'], right_on=['LPERMCO','Year'], how='left')
final_merge = final_merge[['date','PERMNO','PERMCO','EXCHCD','RET','Mkt_Cap','be']]
# drop any na
final_merge = final_merge.dropna(how='any',axis=0) 

# set the year and month 
final_merge['Year'] = final_merge.date//10000
final_merge['Month'] = final_merge['date']//100%100
final_merge = final_merge.loc[final_merge.PERMCO.isin(Non_financial)]

# Adds up stocks from the same company based on PERMCO, 
crsp_mergedBy_JUNE = final_merge[final_merge.Month == 6][['PERMCO', 'Year', 'Mkt_Cap', 'EXCHCD']]

# sorting into deciles
# NYSE stocks are used as breakpoints for all the stocks in NYSE AMEX and NASDAQ
size_bins = crsp_mergedBy_JUNE.loc[crsp_mergedBy_JUNE.EXCHCD==1].groupby('Year')['Mkt_Cap'].apply(lambda x: pd.qcut(x, 10, retbins=True\
                                 ,labels=False, duplicates='drop')[1])  
                
crsp_mergedBy_JUNE['size_decile'] = crsp_mergedBy_JUNE.groupby('Year')['Mkt_Cap']\
        .apply(lambda x: pd.cut(x, bins=size_bins[x.name], labels=False))+1 


crsp_size_merged = pd.merge(final_merge, crsp_mergedBy_JUNE, on=['PERMCO','Year'], how='outer')
# shift by 6 months
crsp_size_merged['size_Rank'] = crsp_size_merged.groupby('PERMNO')['size_decile'].shift(6)
crsp_size_merged['lagged_MktCap'] = crsp_size_merged.groupby('PERMNO')['Mkt_Cap_x'].shift(1)
# remove the NAs and 0s from the lagged market
crsp_size_merged = crsp_size_merged[crsp_size_merged.lagged_MktCap != 0]
crsp_size_merged = crsp_size_merged[(~crsp_size_merged['size_Rank'].isnull()) &
                                    (~crsp_size_merged['lagged_MktCap'].isnull())]

# size return
size_weighted_mean = lambda x: np.average(x, weights=crsp_size_merged.loc[x.index, 'lagged_MktCap'])
size_ret = crsp_size_merged.groupby(['Year','Month','size_Rank'])['RET'].apply(size_weighted_mean).reset_index()





