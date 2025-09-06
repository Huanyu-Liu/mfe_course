n=200;
t=750;

condnum=10;
tau=linspace(1,10,n)';
sigma=diag(tau);

x=randn(t,n);
y=x*sqrtm(sigma);
sample=(y'*y)./t;

sigmabar=sum(diag(sample))/n;
target=sigmabar*eye(n);
denominator=(1/n)*norm(sample-target,'fro')^2;
numerator=(n/t)*sigmabar^2;
a=numerator/denominator;

sigmahat=a*target+(1-a)*sample;

norm(sample-sigma,'fro')^2
norm(sigmahat-sigma,'fro')^2
cond(sample)
cond(sigmahat)

