import numpy as np
import math
class G2Model:
    def __init__(self, phi = 0.03, rho = 0.7, a = 0.1, b = 0.3, sigma = 0.03, eta = 0.08):
        self.phi = phi
        self.rho = rho
        self.a = a
        self.b = b
        self.sigma = sigma
        self.eta = eta

    def bond(self, S, face, x0, y0, path_count):
        delta = 1 / 252
        delta_sqrt = math.sqrt(delta)
        step = int(S * 252)
        x = np.full(path_count,x0)
        y = np.full(path_count,y0)
        R = np.zeros(path_count)
        for i in range(step):
            normal1 = np.random.normal(0, delta_sqrt, path_count)
            normal2 = np.random.normal(0, delta_sqrt, path_count)
            x += -self.a * x * delta + self.sigma * normal2
            y += -self.b * y * delta + self.eta * (self.rho * normal1 + math.sqrt(1 - self.rho * self.rho) * normal2)
            r = x + y + self.phi
            R += r
        return  face * np.exp(-delta * R).mean()

    def option(self,T,strike,S,face,path_count):
        delta = 1 / 252
        delta_sqrt = math.sqrt(delta)
        step = int(T * 252)
        x = np.zeros(path_count)
        y = np.zeros(path_count)
        R = np.zeros(path_count)
        bond = np.zeros(path_count)
        for i in range(step):
            normal1 = np.random.normal(0, delta_sqrt, path_count)
            normal2 = np.random.normal(0, delta_sqrt, path_count)
            x += -self.a * x * delta + self.sigma * normal2
            y += -self.b * y * delta + self.eta * (self.rho * normal1 + math.sqrt(1 - self.rho * self.rho) * normal2)
            r = x + y + self.phi
            R += r
        for i in range(path_count):
            bond[i] = self.bond(S - T, face, x[i], y[i], path_count)
        option = np.maximum(strike - bond, 0) * np.exp(-R * delta)
        return option.mean()