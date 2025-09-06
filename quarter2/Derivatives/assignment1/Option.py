import math
import random
from sklearn.linear_model import LinearRegression
import numpy as np
import pandas as pd

class Option:
    def __init__(self, s0, k, r, sigma, t, is_european, is_call):
        self.s0 = s0
        self.k = k
        self.r = r
        self.sigma = sigma
        self.t = t
        self.is_european = is_european
        self.is_call = is_call
    def set_t(self,t):
        self.t = t
    def binomial_tree(self,h):
        u = math.exp(self.r * h + self.sigma * math.sqrt(h))
        d = math.exp(self.r * h - self.sigma * math.sqrt(h))
        p = (math.exp(self.r * h) - d) / (u - d)
        pd = 1 - p
        stock_price = list()
        payoff = list()
        stock_price.append(self.s0 * math.pow(u, self.t))
        if self.is_call:
            payoff.append(max(stock_price[0] - self.k,0))
        else:
            payoff.append(max(self.k - stock_price[0], 0))
        for i in range(self.t):
            stock_price.append(stock_price[i] * d / u)
            if self.is_call:
                payoff.append(max(stock_price[i + 1] - self.k, 0))
            else:
                payoff.append(max(self.k - stock_price[i + 1], 0))

        if self.is_european:
            for i in range(self.t,0,-1):
                for j in range(i):
                    payoff[j] = payoff[j] * p + payoff[j + 1] * pd


            payoff[0] *= math.exp(-self.r * self.t * h)

        else:
            for i in range(self.t,0,-1):
                for j in range(i):
                    stock_price[j] /= u
                    continue_p = math.exp(-self.r * h) * (p * payoff[j] + pd * payoff[j + 1])
                    if self.is_call:
                        exercise_p = max(stock_price[j] - self.k,0)
                    else:
                        exercise_p = max(self.k - stock_price[j], 0)
                    payoff[j] = max(continue_p,exercise_p)


        return payoff[0]

    def binary(self, h):
        u = math.exp(self.r * h + self.sigma * math.sqrt(h))
        d = math.exp(self.r * h - self.sigma * math.sqrt(h))
        p = (math.exp(self.r * h) - d) / (u - d)
        pd = 1 - p
        stock_price = list()
        payoff = list()
        stock_price.append(self.s0 * math.pow(u, self.t))
        if stock_price[0] > self.k:
            payoff.append(1)
        else:
            payoff.append(0)

        for i in range(self.t):
            stock_price.append(stock_price[i] * d / u)
            if stock_price[i + 1] > self.k:
                payoff.append(1)
            else:
                payoff.append(0)


        for i in range(self.t,0,-1):
            for j in range(i):
                payoff[j] = payoff[j] * p + payoff[j + 1] * pd


        payoff[0] *= math.exp(-self.r)

        return payoff[0]

    def constant_dividend(self,h,dividend,date):
        u = math.exp(self.sigma * math.sqrt(h))
        d = 1 / u
        p = (math.exp(self.r * h) - d) / (u - d)
        pd = 1 - p
        stock_price = list()
        payoff = list()
        dividend_count = len(date)
        stock_price.append(self.s0 * math.pow(u, self.t) * pow(1 - dividend, dividend_count))

        if self.is_call:
            payoff.append(max(stock_price[0] - self.k,0))
        else:
            payoff.append(max(self.k - stock_price[0], 0))
        for i in range(self.t):
            stock_price.append(stock_price[i] * d / u)
            if self.is_call:
                payoff.append(max(stock_price[i + 1] - self.k, 0))
            else:
                payoff.append(max(self.k - stock_price[i + 1], 0))

        if self.is_european:
            for i in range(self.t,0,-1):
                for j in range(i):
                    payoff[j] = payoff[j] * p + payoff[j + 1] * pd


            payoff[0] *= math.exp(-self.r * self.t * h)

        else:
            for i in range(self.t,0,-1):
                for j in range(i):
                    if i in date:
                        stock_price[j] = stock_price[j] / u / (1 - dividend)
                    else:
                        stock_price[j] /= u
                    continue_p = math.exp(-self.r * h) * (p * payoff[j] + pd * payoff[j + 1])
                    if self.is_call:
                        exercise_p = max(stock_price[j] - self.k,0)
                    else:
                        exercise_p = max(self.k - stock_price[j], 0)
                    payoff[j] = max(continue_p,exercise_p)

        return payoff[0]

    def staddle(self,h,dividend,date):
        u = math.exp(self.sigma * math.sqrt(h))
        d = 1 / u
        p = (math.exp(self.r * h) - d) / (u - d)
        pd = 1 - p
        stock_price = list()
        payoff = list()
        dividend_count = len(date)
        stock_price.append(self.s0 * math.pow(u, self.t) * pow(1 - dividend, dividend_count))
        payoff.append(max(stock_price[0] - self.k, self.k - stock_price[0]))
        for i in range(self.t):
            stock_price.append(stock_price[i] * d / u)
            payoff.append(max(stock_price[i + 1] - self.k, self.k - stock_price[i + 1]))

        if self.is_european:
            for i in range(self.t, 0, -1):
                for j in range(i):
                    payoff[j] = payoff[j] * p + payoff[j + 1] * pd

            payoff[0] *= math.exp(-self.r * self.t * h)

        else:
            for i in range(self.t, 0, -1):
                for j in range(i):
                    if i in date:
                        stock_price[j] = stock_price[j] / u / (1 - dividend)
                    else:
                        stock_price[j] /= u
                    continue_p = math.exp(-self.r * h) * (p * payoff[j] + pd * payoff[j + 1])
                    exercise_p = max(stock_price[j] - self.k, self.k - stock_price[j])
                    payoff[j] = max(continue_p, exercise_p)

        return payoff[0]

    def asian(self,h,path_count):
        h_sqrt = math.sqrt(h)

        price = 0
        for i in range(path_count):
            stock_price = self.s0
            sum = self.s0
            for i in range(self.t - 1):
                delta_s = self.r * h * stock_price + self.sigma * stock_price * random.normalvariate(0,h_sqrt)
                stock_price += delta_s
                sum += stock_price

            payoff = max(sum / self.t - self.k,0)
            price += payoff
        return price/path_count

    def lsmc(self,h,path_count):
        discount = math.exp(-self.r * h)
        h_sqrt = math.sqrt(h)
        patha = np.zeros([path_count,self.t])
        tempa = np.zeros(self.t)
        tempa[0] = self.s0
        for i in range(path_count):
            for j in range(self.t - 1):
                delta_sa = tempa[j] * self.r * h + self.sigma * tempa[j] * random.normalvariate(0, h_sqrt)
                tempa[j+1] = tempa[j] + delta_sa
            patha[i] = tempa
        indexa = np.zeros(path_count,dtype=int)
        # path = list()
        #
        # for i in range(path_count):
        #     temp = [self.s0]
        #     for j in range(self.t - 1):
        #         delta_s = temp[j] * self.r * h + self.sigma * temp[j] * random.normalvariate(0,h_sqrt)
        #         temp.append(temp[j] + delta_s)
        #     path.append(temp)
        # index = list()
        # discount = math.exp(-self.r * h)

        for i in range(path_count):
            if patha[i][self.t - 1] < self.k:
                indexa[i] = self.t - 1
            else:
                indexa[i] = self.t
        testa = np.zeros(path_count)
        #indexa = indexa.astype(int)
        for i in range(self.t - 2, 0, -1):
            for j in range(path_count):
                if indexa[j] < self.t:
                    power = self.t - indexa[j]
                    testa[j] = pow(discount,power) * (self.k - patha[j][indexa[j]])

            '''
            matrix calculation
            '''
            df = pd.DataFrame(patha)
            b = np.zeros(3)
            A = np.zeros([3,3])
            b[0] = testa.sum()
            b[1] = (testa * df[i]).sum()
            b[2] = (testa * np.square(df[i])).sum()
            A[0][0] = path_count
            A[0][1] = df[i].sum()
            A[0][2] = np.square(df[i]).sum()
            A[1][1] = A[0][2]
            A[1][2] = (df[i] * np.square(df[i])).sum()
            A[2][2] = np.square(np.square(df[i])).sum()
            A[1][0] = A[0][1]
            A[2][0] = A[0][2]
            A[2][1] = A[1][2]
            a = np.linalg.solve(A,b)
            ecv = a[0] + a[1] * df[i] + a[2] * np.square(df[i])




            # x = pd.concat([df.iloc[:,i],np.square(df.iloc[:,i])],axis=1)
            # a = LinearRegression().fit(x,testa)
            # ecv = a.coef_[0] * x.iloc[:,0] + a.coef_[1] * x.iloc[:,1] + a.intercept_
            for j in range(path_count):
                if (ecv[j] < self.k - patha[j][i]) and (self.k - patha[j][i] > 0):
                    indexa[j] = i
        sum = 0
        for i in range(path_count):
            if indexa[i] < self.t:
                sum += pow(discount,indexa[i]) * (self.k - patha[i][indexa[i]])


        # for i in range(path_count):
        #     if path[i][self.t - 1] < self.k:
        #         index.append(self.t - 1)
        #     else:
        #         index.append(self.t)
        # for i in range(self.t - 2,0,-1):
        #     test = list()
        #     for j in range(path_count):
        #         if index[j] < self.t:
        #             power = self.t - index[j]
        #             test.append(pow(discount,power) * (self.k - path[j][index[j]]))
        #         else:
        #             test.append(0)
        #     df = pd.DataFrame(path)
        #     x = pd.DataFrame({'x':df.iloc[:,i].values,'x2':np.square(df.iloc[:,i].values)})
        #     #print(x.shape)
        #     y = np.array(test)
        #     #print(y.shape)
        #     lm = LinearRegression()
        #     a = lm.fit(x,y)
        #     ecv = a.coef_[0] * df.iloc[:,i] + a.coef_[1] * x.iloc[:,1] + a.intercept_
        #     for j in range(path_count):
        #         if (ecv.values.item(j) < self.k - path[j][i]) and (self.k - path[j][i] > 0):
        #             index[j] = i
        # sum = 0
        #
        # for i in range(path_count):
        #     if index[i] < self.t:
        #         sum += pow(discount,index[i]) * (x - path[i][index[i]])
        return sum / path_count