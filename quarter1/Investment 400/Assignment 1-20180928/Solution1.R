mysoln <- list(student = c("Huanyu Liu", "Meghana Rao", "Chengbo Du", "Shardul Kulkarni"))

#1
invest = 10000
year = 3
ear = 0.06
total_asset_a = invest * (1 + ear) ^ 3

apr_q = 0.06
total_asset_b = invest * (1 + apr_q/4) ^ 12

apr_m = 0.06
total_asset_c = invest * (1 + apr_m/12) ^ 36
total_asset_c

# answers
a = total_asset_a
b = 11956.18
c = 11966.81


#2
payment = 500000 * 0.1
pv = 0
annual_rate = 0.05
for (i in 1:10){
  pv = pv + payment / (1 + annual_rate) ^ (i * 3)
}

pv_b = 700000

function payment(pv, r, n, step = 1){
  a1 = 1 / (1 + r) ^ step
  q = a1
  return pv * (1 - q) / a1 / (1 - q^(n/step))
}

# pv_b = payment / (1 + annual_rate) ^ 3 + payment / (1 + annual_rate) ^ 6 + ... + payment / (1 + annual_rate) ^ 30
# a1 = 1 / (1 + annual_rate) ^ 3, q = 1 / (1 + annual_rate) ^ 3, n = 10
# pv_b = payment * a1 * (1 - q^n) / (1 - q)
# payment = pv_b * (1 - q) / a1 / (1 - q^n)

q = 1 / (1 + annual_rate)^3
a1 = q
payment_b = pv_b * (1 - q) / a1 / (1 - q^10) 

# answer
a = 313279.3
# The present value of this agreement is $313279.3, which is significantly less than $500K. A penny now is more valuable than a penny in the future.
b = 143552.3

#3

year = 35
r = 0.04
fv = 0
save = 200000 * 0.3
for (i in 0:34){
  fv = fv + save * (1 + r) ^ i
}

retirement_year = 20
pv = fv

# pv = consume / (1 + r) + consume / (1 + r)^2 + consume / (1 + r) ^ 3 + ... + consume / (1 + r) ^ 20
# a1 = 1 / (1 + r), q = 1 / (1 + r), n = 20
# consume = pv * (1 - q) / a1 / (1 - q^n)

a1 = 1 / (1 + r)
q = a1
consume = pv * (1 - q) / a1 / (1 - q^retirement_year)

# answer
a = 4419133

b = 325167.6


#4
apr_m = 0.07
ear = (1 + apr_m / 12) ^ 12 - 1

pv = 400000
month = 30 * 12

# pv = payment / (1 + apr_m / 12) + payment / (1 + apr_m / 12)^2 + ... + payment / (1 + apr_m / 12)^360
# a1 = 1 / (1 + apr_m / 12),  