import numpy as np
import math
from scipy.stats import norm

class Option:
    def __init__(self,s0,r,sigma,k,t):
        self.s0 = s0
        self.r = r
        self.sigma = sigma
        self.k = k
        self.t = t

    def set_s0(self,s0):
        self.s0 = s0

    def log_s(self,N,delta_x):
        stock_price = np.zeros(2 * N + 1)
        x0 = math.log(self.s0)
        for i in range(2 * N + 1):
            stock_price[i] = x0 + (N - i) * delta_x
        return stock_price


    def efd(self,N,delta,delta_x):
        temp1 = self.sigma * self.sigma / delta_x / delta_x
        temp2 = (self.r - self.sigma * self.sigma / 2 ) / 2 / delta_x
        m = int(self.t / delta)
        pu = delta * (temp1 / 2 + temp2)
        pm = 1 - delta * temp1 - self.r * delta
        pd = delta * (temp1 / 2 - temp2)
        A = np.full((2 * N + 1, 2 * N + 1), 0.0)
        A[0][0] = pu
        A[0][1] = pm
        A[0][2] = pd
        A[2 * N][2 * N - 2] = pu
        A[2 * N][2 * N - 1] = pm
        A[2 * N][2 * N] = pd
        B = np.zeros(2 * N + 1)
        stock_price = self.log_s(N,delta_x)
        B[2 * N] = math.exp(stock_price[2 * N - 1]) - math.exp(stock_price[2 * N])
        for i in range(1, 2 * N):
            A[i][i - 1] = pu
            A[i][i] = pm
            A[i][i + 1] = pd
        F = np.maximum(self.k - np.exp(stock_price),0)
        for i in range(m):
            F = np.matmul(A,F) + B
        return F

    def ifd(self,N,delta,delta_x):
        temp1 = self.sigma * self.sigma / delta_x / delta_x
        temp2 = (self.r - self.sigma * self.sigma / 2 ) / delta_x
        pu = -0.5 * delta * (temp1 + temp2)
        pm = 1 + delta * temp1 + self.r * delta
        pd = -0.5 * delta * (temp1 - temp2)
        m = int(self.t / delta)
        A = np.full((2 * N + 1, 2 * N + 1), 0.0)
        A[0][0] = 1
        A[0][1] = -1
        A[2 * N][2 * N - 1] = 1
        A[2 * N][2 * N] = -1
        for i in range(1, 2 * N):
            A[i][i - 1] = pu
            A[i][i] = pm
            A[i][i + 1] = pd
        B = np.zeros(2 * N + 1)
        stock_price = self.log_s(N, delta_x)
        B[2 * N] = math.exp(stock_price[2 * N - 1]) - math.exp(stock_price[2 * N])
        F = np.maximum(self.k - np.exp(stock_price), 0)

        for i in range(m):
            B[1:-1] = F[1:-1]
            F = np.linalg.solve(A,B)
        return F

    def cnfd(self,N,delta,delta_x):
        temp1 = self.sigma * self.sigma / delta_x / delta_x
        temp2 = (self.r - self.sigma * self.sigma / 2) / delta_x
        pu = -0.25 * delta * (temp1 + temp2)
        pm = 1 + delta * temp1 / 2 + self.r * delta / 2
        pd = -0.25 * delta * (temp1 - temp2)

        m = int(self.t / delta)
        A = np.full((2 * N + 1, 2 * N + 1), 0.0)
        A[0][0] = 1
        A[0][1] = -1
        A[2 * N][2 * N - 1] = 1
        A[2 * N][2 * N] = -1
        for i in range(1, 2 * N):
            A[i][i - 1] = pu
            A[i][i] = pm
            A[i][i + 1] = pd
        stock_price = self.log_s(N, delta_x)
        F = np.maximum(self.k - np.exp(stock_price), 0)
        B = np.zeros(2 * N + 1)
        B[2 * N] = math.exp(stock_price[2 * N - 1]) - math.exp(stock_price[2 * N])

        for i in range(m):
            B[1:-1] = F[0:-2] * -pu - (pm - 2) * F[1:-1] - pd * F[2:]
            F = np.linalg.solve(A,B)
            #print('B[1]' + str(B[1]) + ' B[-2]' + str(B[-2]))

        return F

    def black_schole(self):
        d1 = 1 / (self.sigma * math.sqrt(self.t)) * (math.log(self.s0 / self.k) + (self.r + self.sigma * self.sigma / 2) * self.t)
        d2 = d1 - self.sigma * math.sqrt(self.t)
        return -self.s0 * norm.cdf(-d1) + self.k * math.exp(- self.r * self.t) * norm.cdf(-d2)