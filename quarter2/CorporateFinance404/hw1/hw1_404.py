import pandas as pd
import numpy as np

path = "/Users/huanyu/Desktop/Corporate Finance 404/hw1"
sp500 = pd.read_csv(path+"/sp500.csv")
daily = pd.read_csv(path + "/daily.csv")
monthly = pd.read_csv(path + "/30days.csv")
DGS1 = pd.read_csv(path + "/DGS1.csv")
DGS5 = pd.read_csv(path + "/DGS5.csv")

def geometric_mean(array,length):
    return np.prod(1 + array) ** (1 / length) - 1

def excess(df,colname):
    sp500_excess = sp500.join(df, how="outer")
    sp500_excess = sp500_excess.dropna()
    sp500_excess['excess'] = sp500_excess['vwretd'] - sp500_excess[colname]
    return sp500_excess

sp500['caldt'] = pd.to_datetime(sp500['caldt'],format="%Y%m%d").dt.date
sp500.set_index('caldt',inplace=True)

daily['qdate'] = pd.to_datetime(daily['qdate'],format="%Y%m%d").dt.date
daily.set_index('qdate',inplace=True)
daily = daily[['ffefrt']] / 360 / 100

monthly['CALDT'] = pd.to_datetime(monthly['CALDT']).dt.date
monthly.set_index('CALDT',inplace=True)
monthly = monthly[['TDBIDYLD']]

DGS1.set_index(pd.to_datetime(DGS1['DATE']).dt.date,inplace=True)
DGS1 = DGS1[['DGS1']]
DGS5.set_index(pd.to_datetime(DGS5['DATE']).dt.date,inplace=True)
DGS5 = DGS5[['DGS5']]
DGS1 = DGS1[DGS1['DGS1'] != '.']
DGS1['DGS1'] = DGS1['DGS1'].astype(float) / 360 / 100
DGS5 = DGS5[DGS5['DGS5'] != '.']
DGS5['DGS5'] = DGS5['DGS5'].astype(float) / 360 / 100

sp500_excess_daily = excess(daily,'ffefrt')
sp500_excess_monthly = excess(monthly,'TDBIDYLD')
sp500_excess_1year = excess(DGS1,'DGS1')
sp500_excess_5year = excess(DGS5, 'DGS5')

sp500_arith_daily = np.mean(sp500['vwretd'])
sp500_geo_daily = geometric_mean(sp500['vwretd'],len(sp500['vwretd']))

excess_arith_daily = np.mean(sp500_excess_daily['excess'])
excess_geo_daily = geometric_mean(sp500_excess_daily['excess'],len(sp500_excess_daily['excess']))

s = pd.Series(np.arange(sp500.shape[0]), index=pd.to_datetime(sp500.index))
month_end = s.resample("M").max().values
month_end = np.insert(month_end, 0,-1)
year_end = s.resample("Y").max().values
year_end = np.insert(year_end, 0 , -1)

sp500_monthly = list()

for i in range(month_end.size - 1):
    sp500_monthly.append(np.prod(sp500.iloc[month_end[i] + 1: month_end[i+1]] + 1) - 1)

sp500_monthly = np.array(sp500_monthly)
sp500_arith_month = np.mean(sp500_monthly)
sp500_geo_month = geometric_mean(sp500_monthly, sp500_monthly.size)

excess_monthly = list()

for i in range(month_end.size - 1):
    excess_monthly.append(np.prod(sp500_excess_monthly.iloc[month_end[i] + 1: month_end[i+1],2] + 1) - 1)

excess_monthly = np.array(excess_monthly)
excess_arith_monthly = np.mean(excess_monthly)
excess_geo_monthly = geometric_mean(excess_monthly, excess_monthly.size)

sp500_annually = list()
for i in range(year_end.size - 1):
    sp500_annually.append(np.prod(sp500.iloc[year_end[i] + 1: year_end[i+1]] + 1) - 1)
sp500_annually = np.array(sp500_annually)
sp500_arith_annual = np.mean(sp500_annually)
sp500_geo_annual = geometric_mean(sp500_annually,sp500_annually.size)

excess_annually = list()
for i in range(year_end.size - 1):
    excess_annually.append(np.prod(sp500_excess_1year.iloc[year_end[i] + 1: year_end[i+1],2] + 1) - 1)
excess_annually = np.array(excess_annually)

excess_arith_annually = np.mean(excess_annually)
excess_geo_annually = geometric_mean(excess_annually,excess_annually.size)

sp500_5year = list()

for i in range(0,year_end.size -3, 5):
    sp500_5year.append(np.prod(sp500.iloc[year_end[i] + 1: year_end[i+5]] + 1) - 1)

sp500_5year = np.array(sp500_5year)
sp500_arith_5year = np.mean(sp500_5year)
sp500_geo_5year = geometric_mean(sp500_5year,sp500_5year.size)

excess_5year = list()

for i in range(0,year_end.size - 3,5):
    excess_5year.append(np.prod(sp500_excess_5year.iloc[year_end[i] + 1: year_end[i+5],2] + 1) - 1)
excess_5year = np.array(excess_5year)

excess_arith_5year = np.mean(excess_5year)
excess_geo_5year = geometric_mean(excess_5year, excess_5year.size)