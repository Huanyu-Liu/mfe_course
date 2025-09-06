from bond import Bond
from g2_model import G2Model
import time

start = time.time()
r0 = 0.05
sigma = 0.18
k = 0.82
r_bar = 0.05
face = 1000
T = 0.5
path_count = 1000
bond = Bond(r0,sigma,k,r_bar)
coupon = list()
coupon_T = list()
for i in range(7):
    coupon.append(30)
    coupon_T.append((i + 1) * 0.5)
coupon.append(1030)
coupon_T.append(4)
print('Q1(a): {}'.format(bond.price(T,face,path_count)))
print('Q1(b): {}'.format(bond.coupon_bond(coupon,coupon_T,path_count)))

strike = 980
T = 0.25
S = 0.5
print('Q1(c): {}'.format(bond.zcb_option(T,strike,S,face,path_count)))
print('Q1(d): {}'.format(bond.cb_option(T,strike,coupon_T,coupon,path_count)))

T = 0.5
S = 1
k = 0.92
r_bar = 0.055
bond2 = Bond(r0, sigma, k, r_bar)
print('Q2(a): {}'.format(bond2.cir_option(T,strike,S,face,path_count)))
#print('explicit')
print('Q2(b): {}'.format(bond2.explicit_cir_option(T, 0.98, S, face, path_count)))

strike = 985
bond3 = G2Model()
print('Q3: {}'.format(bond3.option(T, strike, S, face, path_count)))

end = time.time()
print('Total runtime: {}'.format(end - start))