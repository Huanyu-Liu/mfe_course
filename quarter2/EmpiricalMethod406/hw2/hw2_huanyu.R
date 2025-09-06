r = vector()
r[1] = 0
r[2] = 1
phi0 = 2
phi1 = 1.1
phi2 = -0.25
sigma = 1
for (i in c(2:1000)){
  r[i+1] = phi0 + phi1 * r[i] + phi2 * r[i-1] + rnorm(1,0,sigma)
}
acf(r, xlim = c(0,20))
