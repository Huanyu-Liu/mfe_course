#%%
import numpy as np
import scipy.linalg as lg
import pandas as pd
#%%
index = ['AAA', 'AA', 'A','BBB', 'BB', 'B', 'CCC','Default']
p0 = np.zeros(shape=(8,8))
for i in range(8):
    p0[i][i] = 1
p0_df = pd.DataFrame(p0,index=index, columns=index)
p1 = np.array(
    [[90.81,8.33,0.68,0.06,0.12,0,0,0],
     [0.7,90.65,7.79,0.64,0.06,0.14,0.02,0],
     [0.09,2.27,91.05,5.52,0.74,0.26,0.01,0.06],
     [0.02,0.33,5.95,86.93,5.3,1.17,1.12,0.18],
     [0.03,0.14,0.67,7.73,80.53,8.84,1,1.06],
     [0,0.11,0.24,0.43,6.48,83.46,4.07,5.2],
     [0.22,0,0.22,1.3,2.38,11.24,64.86,19.79],
     [0,0,0,0,0,0,0,100]]
) / 100
p1_df = pd.DataFrame(p1,index=index,columns=index)
print('P0')
print(p0_df)
print('\n\nP1')
print(p1_df)

lamb = lg.logm(p1)
lambda_df = pd.DataFrame(lamb,index=index, columns=index)
print(lambda_df)
#%%
default = np.zeros(shape=(7,7))
for i,v in enumerate([1,2,3,4,5,7,10]):
    p = lg.expm(lamb * v)
    default[i] = p[:7, 7]
default_df = pd.DataFrame(default.T,index=index[:-1],columns=[1,2,3,4,5,7,10])
print(default_df)
#%%
bond_price = 0
for i in range(12):
    default_rate = lg.expm(lamb * (i + 1) * 0.5)[3][-1]
    bond_price += 3 * (1 - default_rate)
bond_price += 100 * (1 - default_rate) + 60 * default_rate
print(bond_price)
#%%
recovery = 0.6
default3 = lg.expm(lamb * 3)[3][-1]
default5 = lg.expm(lamb * 5)[3][-1]
default10 = lg.expm(lamb * 10)[3][-1]
cds_spread3 = (1 - recovery) * default3 / (1 - default3)
cds_spread5 = (1 - recovery) * default5 / (1 - default5)
cds_spread10 = (1 - recovery) * default10 / (1 - default10)
print(cds_spread3)
print(cds_spread5)
print(cds_spread10)