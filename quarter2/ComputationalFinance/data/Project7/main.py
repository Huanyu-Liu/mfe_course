import math
from Option import Option
from American_option import American_option
import numpy as np
import pandas as pd

path = '/Users/huanyu/Desktop/ComputationalFinance/data/Project7/'

# Problem 1
s0 = 10
sigma = 0.2
r = 0.04
x0 = math.log(s0)
delta = 0.002
delta_x = sigma * math.sqrt(delta)
k = 10
T = 0.5
N = 200
op = Option(s0,r,sigma,k,T)
stock_price = np.zeros(13)
efd = np.zeros(13)
ifd = np.zeros(13)
cnfd = np.zeros(13)
black_scholes = np.zeros(13)
print("Problem 1: \n")
for i in (1,3,4):
    delta_x = sigma * math.sqrt(i * delta)
    price_efd = op.efd(N, delta, delta_x)
    price_ifd = op.ifd(N, delta, delta_x)
    price_cnfd = op.cnfd(N, delta, delta_x)

    stock_price = np.zeros(13)
    efd = np.zeros(13)
    ifd = np.zeros(13)
    cnfd = np.zeros(13)
    black_scholes = np.zeros(13)
    for j in range(4, 17):
        temp = int(round((math.log(j) - x0) / delta_x))
        op.set_s0(j)
        stock_price[j - 4] = j
        efd[j - 4] = price_efd[N - temp]
        ifd[j - 4] = price_ifd[N - temp]
        cnfd[j - 4] = price_cnfd[N - temp]
        bs_price = op.black_schole()
        black_scholes[j - 4] = bs_price

    df = pd.DataFrame({'stock price':stock_price,'efd':efd,'ifd':ifd,'cnfd':cnfd,'black scholes':black_scholes})
    df['efd_error'] = df['efd'] / df['black scholes'] - 1
    df['ifd_error'] = df['ifd'] / df['black scholes'] - 1
    df['cnfd_error'] = df['cnfd'] / df['black scholes'] - 1
    df.to_csv(path + 'problem1_' + str(i) + '.csv',index=False)
    print('\nDelta_x = sigma * sqrt(' + str(i) + ' * delta_t):\n')
    print(df)

    op.set_s0(s0)


# Problem 2
is_call = True
delta_s = 0.5
N = int(round(s0 / delta_s))

american_call = American_option(s0,r,sigma,k,T,is_call)
american_put = American_option(s0,r,sigma,k,T,not is_call)
american_call_efd = american_call.efd(delta,delta_s)
american_put_efd = american_put.efd(delta,delta_s)
print('\nProblem 2\n')

call_efd = np.zeros(13)
call_ifd = np.zeros(13)
call_cnfd = np.zeros(13)
put_efd = np.zeros(13)
put_ifd = np.zeros(13)
put_cnfd = np.zeros(13)

for i in (1,4,5):
    delta_s = 0.25 * i
    N = int(round(s0 / delta_s))
    american_call_efd = american_call.efd(delta, delta_s)
    american_put_efd = american_put.efd(delta, delta_s)
    american_call_ifd = american_call.ifd(delta, delta_s)
    american_put_ifd = american_put.ifd(delta, delta_s)
    american_call_cnfd = american_call.cnfd(delta, delta_s)
    american_put_cnfd = american_put.cnfd(delta, delta_s)


    print('\nDelta S = ' + str(i * 0.25) + ':\n')
    for j in range(4, 17):
        temp = int(round((j - s0) / delta_s))
        call_efd[j - 4] = american_call_efd[N - temp]
        call_ifd[j - 4] = american_call_ifd[N - temp]
        call_cnfd[j - 4] = american_call_cnfd[N - temp]
        put_efd[j - 4] = american_put_efd[N - temp]
        put_ifd[j - 4] = american_put_ifd[N - temp]
        put_cnfd[j - 4] = american_put_cnfd[N - temp]

    df_call = pd.DataFrame({'stock price':stock_price,'efd':call_efd,'ifd':call_ifd,'cnfd':call_cnfd})
    df_put = pd.DataFrame({'stock price':stock_price,'efd':put_efd,'ifd':put_ifd,'cnfd':put_cnfd})
    df_call.to_csv(path + 'problem2call_' + str(i) + '.csv',index=False)
    df_put.to_csv(path + 'problem2put_' + str(i) + '.csv',index=False)
    print('Call:\n')
    print(df_call)
    print('\nPut:\n')
    print(df_put)
    