import numpy as np
import math
from scipy.stats import ncx2
class Bond:
    def __init__(self, r0, sigma, k, r_bar):
        self.r0 = r0
        self.sigma = sigma
        self.k = k
        self.r_bar = r_bar

    def price(self, T, face, path_count):
        delta = 1 / 252
        delta_sqrt = math.sqrt(delta)
        step = int(T * 252)
        rate = np.full(path_count, self.r0)
        R = np.zeros(path_count)
        for i in range(step):
            normal = np.random.normal(0,delta_sqrt,path_count)
            rate += self.k * (self.r_bar - rate) * delta + self.sigma * normal
            R += rate
        return face * np.exp(-R * delta).mean()

    def coupon_bond(self, coupon, T, path_count):
        coupon_count = len(coupon)
        sum = np.zeros(path_count)
        step = int(T[-1] * 252)
        delta = 1 / 252
        delta_sqrt = math.sqrt(delta)
        coupon_step = [int(x * 252) - 1 for x in T]
        rate = np.full(path_count, self.r0)
        R = np.zeros(path_count)
        counter = 0
        for i in range(step):
            normal = np.random.normal(0, delta_sqrt, path_count)
            rate += self.k * (self.r_bar - rate) * delta + self.sigma * normal
            R += rate
            if i in coupon_step:
                sum += coupon[counter] * np.exp(-R * delta)
                counter += 1
        return sum.mean()

    def zcb_option(self,T,strike,S,face,path_count):
        delta = 1 / 252
        delta_sqrt = math.sqrt(delta)
        step = int(T * 252)
        rate = np.full(path_count, self.r0)
        R_option = np.zeros(path_count)
        B = 1 / self.k * (1 - math.exp(-self.k * (S - T)))
        A = math.exp((self.r_bar - self.sigma * self.sigma / 2 / self.k / self.k) * (B - S + T) - self.sigma * self.sigma / 4 / self.k * B * B)
        for i in range(step):
            normal = np.random.normal(0,delta_sqrt,path_count)
            rate += self.k * (self.r_bar - rate) * delta + self.sigma * normal
            R_option += rate

        bond = face * A * np.exp(-B * rate)
        option = np.maximum(bond - strike, 0) * np.exp(-R_option * delta)
        return option.mean()

    def cb_option(self, T, strike, S, coupon, path_count):
        step = int(T * 252)
        delta = 1 / 252
        original_r = self.r0
        delta_sqrt = math.sqrt(delta)
        coupon_step = [int(x * 252) - 1 for x in S]
        rate = np.full(path_count, self.r0)
        R = np.zeros(path_count)
        S = [x - T for x in S]
        bond = np.zeros(path_count)
        counter = 0
        for i in range(step):
            normal = np.random.normal(0, delta_sqrt, path_count)
            rate += self.k * (self.r_bar - rate) * delta + self.sigma * normal
            R += rate
        for i in range(path_count):
            self.r0 = rate[i]
            bond[i] = self.coupon_bond(coupon,S, path_count)
        option = np.maximum(bond - strike, 0) * np.exp(-R * delta)
        self.r0 = original_r
        return option.mean()

    def cir_bond(self, rate, T, face, path_count):
        delta = 1 / 252
        delta_sqrt = math.sqrt(delta)
        step = int(T * 252)
        bond = np.zeros(path_count)
        for i in range(path_count):
            rate_array = np.full(path_count, rate[i])
            R = np.zeros(path_count)
            for j in range(step):
                normal = np.random.normal(0, delta_sqrt, path_count)
                rate_array += self.k * (self.r_bar - rate_array) * delta + np.sqrt(abs(rate_array)) * self.sigma * normal
                R += rate_array
            bond[i] = face * np.exp(-R * delta).mean()
        return bond


    def cir_option(self,T,strike,S,face,path_count):
        step = int(T * 252)
        delta = 1 / 252
        original_r = self.r0
        delta_sqrt = math.sqrt(delta)
        rate = np.full(path_count, self.r0)
        R = np.zeros(path_count)
        bond = np.zeros(path_count)
        counter = 0
        for i in range(step):
            normal = np.random.normal(0, delta_sqrt, path_count)
            rate += self.k * (self.r_bar - rate) * delta + np.sqrt(rate) * self.sigma * normal
            R += rate

        bond = self.cir_bond(rate, S - T, face, path_count)
        #print(bond)
        option = np.maximum(bond - strike, 0) * np.exp(-R * delta)
        return option.mean()

    def explicit_cir_option(self,T,strike,S,face,path_count):
        theta = math.sqrt(self.k * self.k + 2 * self.sigma * self.sigma)
        phi = 2 * theta / (self.sigma * self.sigma * (math.exp(theta * T) - 1))
        sai = (self.k + theta) / self.sigma / self.sigma
        A, B = self.A_B(T,S)

        rt = self.r0
        bond1 = self.P(0,S, self.r0)
        bond2 = self.P(0,T, self.r0)

        r_star = math.log(A / strike) / B
        x1 = 2 * r_star * (phi + sai + B)
        p1 = 4 * self.k * self.r_bar / self.sigma / self.sigma
        temp = 2 * phi * phi * self.r0 * math.exp(theta * T)
        q1 = temp / (phi + sai + B)
        x2 = 2 * r_star * (phi + sai)
        q2 = temp / (phi + sai)
        chi_1 = ncx2.cdf(x1, p1, q1)
        chi_2 = ncx2.cdf(x2, p1, q2)
        return face * (bond1 * chi_1 - strike * bond2 * chi_2)


    def A_B(self, t, T):
        h1 = math.sqrt(self.k * self.k + 2 * self.sigma * self.sigma)
        h2 = (self.k + h1) / 2
        h3 = 2 * self.k * self.r_bar / self.sigma / self.sigma
        a = h1 * math.exp(h2 * (T - t))
        b = math.exp(h1 * (T - t)) - 1
        c = h2 * b + h1
        return math.pow(a / c, h3), b / c

    def P(self,t, T, rt):
        A, B = self.A_B(t, T)

        return A * math.exp(-B * rt)