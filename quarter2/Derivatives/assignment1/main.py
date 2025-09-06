from Option import Option
import matplotlib.pyplot as plt

s0 = 100
k = 90
r = 0.02
t = 4
h = 0.25
sigma = 0.2
is_call = True
is_european = True

op_call1 = Option(s0,s0,r,sigma,t,is_call,is_european)
op_put1 = Option(s0,s0,r,sigma,t,not is_call, is_european)
# Problem 1
# (a)

straddle1 = op_call1.binomial_tree(h) + op_put1.binomial_tree(h)
print(straddle1)

# (b)

t = 40
h = 0.025
op_call2 = Option(s0,s0,r,sigma,t,is_call,is_european)
op_put2 = Option(s0,s0,r,sigma,t,not is_call,is_european)

straddle2 = op_call2.binomial_tree(h) + op_put2.binomial_tree(h)
print(straddle2)

# (c)

t = 4
h = 0.25

op_binary = Option(s0,k,r,sigma,t,is_call,is_european)

binary = op_binary.binary(h)
print(binary)

# Problem 2
s0 = 10
r = 0.01
t = 250
k = 10
sigma = 0.15
american_call = Option(s0,k,r,sigma,t,is_call, not is_european).binomial_tree(h)
american_put = Option(s0,k,r,sigma,t,not is_call, not is_european).binomial_tree(h)

print(american_call)
print(american_put)

# Problem 3
# (a)

s0 = 10
k = 10
r = 0.02
h = 1/365
t = 200
date = [50,100,150]
dividend = 0.05
sigma = 0.2
op3a1 = Option(s0,k,r,sigma,t,not is_european,is_call)
op3a2 = Option(s0,k,r,sigma,t,not is_european, not is_call)

american_call_d = op3a1.constant_dividend(h,dividend, date)
american_put_d = op3a2.constant_dividend(h,dividend,date)
american_straddle = op3a1.staddle(h,dividend,date)

print(american_call_d)
print(american_put_d)
print(american_straddle)

# Problem 4
# (a)

s0 = 200
r = 0.02
sigma = 0.2
k = 220
t = 365
h = 1/365
path_count = 100000

op4a = Option(s0,k,r,sigma,t,is_european,is_call)
asian_op = op4a.asian(h,path_count)
print(asian_op)

# Problem 5
# (a)

s0 = 200
r = 0.1
sigma = 0.3
k = 220
t = 250
h = 1/t

op5a = Option(s0,k,r,sigma,t,is_european,is_call)
american_put_lsmc = op5a.lsmc(h,path_count)
print(american_put_lsmc)

# (b)

paths = [10,100,1000,10000,100000]
price = list()
for i in paths:
    price.append(op5a.lsmc(h,i))
plt.plot(paths,price)

Ns = [3,10,100,250,1000]
price = list()
for i in Ns:
    op5a.set_t(i)
    price.append(op5a.lsmc(h,path_count))

plt.plot(Ns,price)
plt.show()