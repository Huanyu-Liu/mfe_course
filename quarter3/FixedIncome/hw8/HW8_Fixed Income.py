# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 20:15:43 2019

@author: Serena Peng
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st



if __name__ == '__main__':
    
    #Import the data from last homework 
    path_hw7 = 'S:/UCLA/408 Fixed Income/HW7_String Model/'
    
    #Corrleation
    corr = pd.read_csv(path_hw7 + 'Homework 7 corrin.csv', header = None)

    #Cholesky Decomp of Correlation: Corr = chol@chol.T
    chol = pd.read_csv(path_hw7 + 'Homework 7 corchol.csv', header = None)
    
    #Discount rate: D(t)
    discount = pd.read_excel(path_hw7 + 'Homework 7 pfilea.xlsx', header = None)\
                    .rename(columns = {0: 'Dt'})
    
    #Volatility
    vol = pd.read_excel(path_hw7 + 'Homework 7 sigma.xlsx', header = None)\
                    .to_numpy().reshape((len(discount) - 1))


    '''Question 1: Caps'''
    t = [2, 3, 4, 5, 7, 10] #maturity
    
    #Calculate the strike price
    cms = [] #CMS rate, aka, par rate, aka the strike price
    for i in range(0, len(t)):
        tt = t[i]
        cms.append(2*(1 - discount.loc[tt/0.5-1, 'Dt'])/\
                   np.sum(discount.loc[:tt/0.5-1, 'Dt']))
    
    
#    '''Method 1: Black model'''    
#    #Calculate the forward rate now
#    fwds = ((discount['Dt'].shift()/discount['Dt']-1)*2)
#    
#    caps_black = []
#    for i in range(0, len(t)):
#        k = cms[i]
#        tt = np.arange(0.5, t[i], 0.5) #ignore the first caplet
#        
#        #black model
#        d = (np.log(fwds[1:t[i]*2]/k) + (vol[:t[i]*2-1]**2)*tt/2)/(vol[:t[i]*2-1]*np.sqrt(tt))
#        nd = st.norm.cdf(d)
#        nd_ = st.norm.cdf(d - vol[:t[i]*2-1]*np.sqrt(tt))
#        
#        caplets = (0.5*365/360)*discount.iloc[1:t[i]*2,0]*(fwds[1:t[i]*2]*nd - k*nd_)
#        caps_black.append(np.sum(caplets))
        
    
    '''Method 2: Simulation using String Model (from Homework 7)'''
    #Simulation
    dt = 0.5
    m = 10000
    T = 20 #10-years (semi-annual)
    ft = int(t[-1]/dt)

    #Model the discount rate
    caps_df = np.zeros((m, len(t)))
    for mm in range(0, m):
    
        np.random.seed(mm**2)
        z = (np.random.normal(size = int(ft*T)).reshape((T, ft)).T@chol.to_numpy()[:, :T]).T
        
        dtpath = np.zeros((T, ft+1))
        dtpath[:, 0] =  discount.iloc[:T, 0]
        
        #Construct the whole panel
        for c in range(1, ft+1):
            r = ((1/dtpath[c-1, c-1]) - 1)*2
            
            dtpath[c-1, c] = 1
            dtpath[c:, c] = dtpath[c:, c-1] + r*dtpath[c:, c-1]*dt + \
                    vol[:T-c]*dtpath[c:, c-1]*np.sqrt(dt)*z[c:, c-1]
       
        #Extract only the 6-month discount rate D(t, t+0.5) for each period ahead
        dt05 = np.diag(dtpath)
        ll = (1/dt05[1:] - 1)*2 #drop the 1st rate
        
        #Calculate the payoffs of all caplets for each cap with different strike price
        for i in range(0, len(t)):
            k = cms[i]
            tt = t[i]
            payoff = (0.5*365/360)*np.maximum(ll[:tt*2-1] - k, 0)\
                            *np.cumprod(dt05[:tt*2])[1:]  
            caps_df[mm, i] = np.sum(payoff)
                
    caps_sim = np.mean(caps_df, axis = 0)
    
    output1 = pd.DataFrame({'Strike': cms, 'Caps': caps_sim}, index = t)
    print(output1)
    print('----------------------------------------------------------------')
    
    '''Question 2: at-the-money-forward Swaption
    Step 1: Using simulations to get the forward par rates (strike price)
    Step 2: Get the forward par rate for each path
    Step 3: For each path, compare the par rate with strike, if higher, exercise the swaption
            payoff = fixed leg value with strike par rate - 1
    Step 4: discount the payoff back to time 0'''
    dt = 0.5
    m = 10000
    T = 20 #10-years (semi-annual)
    ft = int(t[-1]/dt)
    
    '''Step 1 & 2'''
    dt05 = np.zeros((m, T))
    
    #Model the discount rate --> find the strike price (forward par rate)
    par1_t = [1, 2, 3, 4] #the duration of swap
    par1_df = np.zeros((m, len(par1_t)))
    fwd1_dt = np.zeros((m, 2*par1_t[-1])) #record the forward d(t)
    
    par2_t = [1, 2, 3]
    par2_df = np.zeros((m, len(par2_t)))
    fwd2_dt = np.zeros((m, 2*par2_t[-1])) #record the forward d(t)

    
    par5_t = [1, 2, 5]
    par5_df = np.zeros((m, len(par5_t)))
    fwd5_dt = np.zeros((m, 2*par5_t[-1])) #record the forward d(t)
 
    for mm in range(0, m):
    
#        np.random.seed((mm*2)**2)
        np.random.seed((mm*2)**2 + 1)
        z = (np.random.normal(size = int(ft*T)).reshape((T, ft)).T@chol.to_numpy()[:, :T]).T
        
        dtpath = np.zeros((T, ft+1))
        dtpath[:, 0] =  discount.iloc[:T, 0]
        
        #Construct the whole panel
        for c in range(1, ft+1):
            r = ((1/dtpath[c-1, c-1]) - 1)*2
            
            dtpath[c-1, c] = 1
            dtpath[c:, c] = dtpath[c:, c-1] + r*dtpath[c:, c-1]*dt + \
                    vol[:T-c]*dtpath[c:, c-1]*np.sqrt(dt)*z[c:, c-1]
    
        #Document the forward D(0.5)
        dt05[mm, :] = np.diag(dtpath)
        
        #solve 1-year forward 1,2,5-year par rate
        par1 = [] #par rate
        for p in par1_t:
            par1.append(2*(1 - dtpath[(1+p)*2-1, 2])/np.sum(dtpath[2:(1+p)*2,2]))
        par1_df[mm, :] = par1
        fwd1_dt[mm, :] = dtpath[2:(1+p)*2,2]
        
        #solve 2-year forward 1,2,5-year par rate
        par2 = [] #par rate
        for p in par2_t:
            par2.append(2*(1 - dtpath[(2+p)*2-1, 4])/np.sum(dtpath[4:(2+p)*2,4]))
        par2_df[mm, :] = par2
        fwd2_dt[mm, :] = dtpath[4:(2+p)*2,4]
    
        #solve 5-year forward 1,2,5-year par rate
        par5 = [] #par rate
        for p in par5_t:
            par5.append(2*(1 - dtpath[(5+p)*2-1,10])/np.sum(dtpath[10:(5+p)*2,10]))
        par5_df[mm, :] = par5
        fwd5_dt[mm, :] = dtpath[10:(5+p)*2,10]
        
    k1 = np.mean(par1_df, axis = 0) #strike price
    k2 = np.mean(par2_df, axis = 0) #strike price
    k5 = np.mean(par5_df, axis = 0) #strike price
    
    par_df = np.concatenate((par1_df, par2_df, par5_df), axis = 1)
    
    strikes =  np.append(np.append(k1, k2), k5)
    
    '''Step 3'''
    fwd_t = [1]*len(par1_t) + [2]*len(par2_t) + [5]*len(par5_t)
    s_t = par1_t + par2_t + par5_t
    
    combo = [(fwd_t[i], s_t[i]) for i in range(0, len(s_t))]
    
    swaption_df = np.zeros((m, len(s_t)))
    for c in range(0, len(combo)):
            
        #strike par rate
        kk = strikes[c] 
        
        #Info about the ff_to_ss swaption
        ff, ss = combo[c]
        
        #all paths of par rate
        pars = par_df[:, c]
        
        #Cashflows
        cf = np.zeros(ss*2) + kk/2
        cf[-1] += 1
        
        if ff == 1:
            payoff = np.sum(fwd1_dt[pars < kk, :(ss*2)]*cf, axis = 1) - 1
            swaption_df[pars < kk, c] = payoff*np.prod(dt05[pars < kk, :(ff*2)], axis = 1)
        elif ff == 2: 
            payoff = np.sum(fwd2_dt[pars < kk, :(ss*2)]*cf, axis = 1) - 1
            swaption_df[pars < kk, c] = payoff*np.prod(dt05[pars < kk, :(ff*2)], axis = 1)
        elif ff == 5: 
            payoff = np.sum(fwd5_dt[pars < kk, :(ss*2)]*cf, axis = 1) - 1
            swaption_df[pars < kk, c] = payoff*np.prod(dt05[pars < kk, :(ff*2)], axis = 1)
    
    #Swaption is the average payoff discounted to now
    swaptions = np.mean(swaption_df, axis = 0)
    
    #Summary
    output2 = pd.DataFrame({'Strike': strikes, 'Swaption': swaptions}, 
                            index = ['1_into_1', '1_into_2', '1_into_3', '1_into_4',
                                     '2_into_1', '2_into_2', '2_into_3', 
                                     '5_into_1', '5_into_2', '5_into_5'])
    print(output2)

    
    