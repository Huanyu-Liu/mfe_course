import numpy as np
import math

class American_option:
    def __init__(self,s0,r,sigma,k,t,is_call = True):
        self.s0 = s0
        self.r = r
        self.sigma = sigma
        self.k = k
        self.t = t
        self.is_call = is_call

    def stock(self,delta_s):
        N = int(self.s0 / delta_s)
        stock_price = np.zeros(2 * N + 1)
        for i in range(2 * N + 1):
            stock_price[i] = self.s0 + (N - i) * delta_s
        return stock_price

    def efd(self,delta,delta_s):
        N = int(self.s0 / delta_s)
        A = np.full((2 * N + 1, 2 * N + 1), 0.0)
        B = np.zeros(2 * N + 1)
        stock_price = self.stock(delta_s)
        m = int(self.t / delta)
        if self.is_call:
            B[0] = stock_price[0] - stock_price[1]
            exercise = np.maximum(stock_price - self.k, 0)
            F = np.copy(exercise)
        else:
            B[2 * N] = stock_price[2 * N - 1] - stock_price[2 * N]
            exercise = np.maximum(self.k - stock_price, 0)
            F = np.copy(exercise)
        pu = np.zeros(2 * N - 1)
        pm = np.zeros(2 * N - 1)
        pd = np.zeros(2 * N - 1)
        for i in range(1,2 * N):
            temp1 = self.r * i / 2
            temp2 = self.sigma * self.sigma * i * i
            pu[i - 1] = delta * (temp1 + temp2 / 2)
            pm[i - 1] = 1 - delta * (temp2 + self.r)
            pd[i - 1] = delta * (-temp1 + temp2 / 2)

        A[0][0] = pu[-1]
        A[0][1] = pm[-1]
        A[0][2] = pd[-1]
        A[2 * N][2 * N - 2] = pu[0]
        A[2 * N][2 * N - 1] = pm[0]
        A[2 * N][2 * N] = pd[0]

        for i in range(1, 2 * N):
            A[i][i - 1] = pu[2 * N - 1 - i]
            A[i][i] = pm[2 * N - 1 - i]
            A[i][i + 1] = pd[2 * N - 1 - i]


        for i in range(m):
            ecv = np.matmul(A,F) + B
            F = np.maximum(exercise,ecv)
        return F

    def ifd(self,delta,delta_s):
        N = int(self.s0 / delta_s)
        A = np.full((2 * N + 1, 2 * N + 1), 0.0)
        B = np.zeros(2 * N + 1)
        stock_price = self.stock(delta_s)
        m = int(self.t / delta)
        if self.is_call:
            B[0] = stock_price[0] - stock_price[1]
            exercise = np.maximum(stock_price - self.k, 0)
            F = np.copy(exercise)
        else:
            B[2 * N] = stock_price[2 * N - 1] - stock_price[2 * N]
            exercise = np.maximum(self.k - stock_price, 0)
            F = np.copy(exercise)
        a1 = np.zeros(2 * N - 1)
        a2 = np.zeros(2 * N - 1)
        a3 = np.zeros(2 * N - 1)
        for i in range(1,2 * N):
            temp1 = self.r * i / 2
            temp2 = self.sigma * self.sigma * i * i
            a1[i - 1] = (temp1 + temp2 / 2) * -delta
            a2[i - 1] = 1 + delta * (temp2 + self.r)
            a3[i - 1] = -delta * (temp2 / 2 - temp1)

        A[0][0] = 1
        A[0][1] = -1
        A[2 * N][2 * N - 1] = 1
        A[2 * N][2 * N] = -1

        for i in range(1, 2 * N):
            A[i][i - 1] = a1[2 * N - 1 - i]
            A[i][i] = a2[2 * N - 1 - i]
            A[i][i + 1] = a3[2 * N - 1 - i]


        for i in range(m):
            B[1:-1] = F[1:-1]
            ecv = np.linalg.solve(A,B)
            F = np.maximum(exercise,ecv)
        return F

    def cnfd(self,delta,delta_s):
        N = int(self.s0 / delta_s)
        A = np.full((2 * N + 1, 2 * N + 1), 0.0)
        B = np.zeros(2 * N + 1)
        stock_price = self.stock(delta_s)
        m = int(self.t / delta)
        if self.is_call:
            B[0] = stock_price[0] - stock_price[1]
            exercise = np.maximum(stock_price - self.k, 0)
            F = np.copy(exercise)
        else:
            B[2 * N] = stock_price[2 * N - 1] - stock_price[2 * N]
            exercise = np.maximum(self.k - stock_price, 0)
            F = np.copy(exercise)
        a1 = np.zeros(2 * N - 1)
        a2 = np.zeros(2 * N - 1)
        a3 = np.zeros(2 * N - 1)
        b1 = np.zeros(2 * N - 1)
        b2 = np.zeros(2 * N - 1)
        b3 = np.zeros(2 * N - 1)
        for i in range(1,2 * N):
            temp1 = self.r * i
            temp2 = self.sigma * self.sigma * i * i
            a1[i - 1] = (temp1 + temp2) / 4 * -delta
            a2[i - 1] = 1 + delta * (temp2 + self.r) / 2
            a3[i - 1] = -delta * (temp2 - temp1) / 4
            b1[i - 1] = -a1[i - 1]
            b2[i - 1] = 2 - a2[i - 1]
            b3[i - 1] = -a3[i - 1]

        A[0][0] = 1
        A[0][1] = -1
        A[2 * N][2 * N - 1] = 1
        A[2 * N][2 * N] = -1

        for i in range(1, 2 * N):
            A[i][i - 1] = a1[2 * N - 1 - i]
            A[i][i] = a2[2 * N - 1 - i]
            A[i][i + 1] = a3[2 * N - 1 - i]


        for i in range(m):
            B[1:-1] = F[0:-2] * b1 + F[1:-1] * b2 + F[2:] * b3
            ecv = np.linalg.solve(A,B)
            F = np.maximum(exercise,ecv)
        return F