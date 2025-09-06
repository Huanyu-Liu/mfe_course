%% Problem 1
% dynamics of St = sigma*St*dWt (no drift term as r = 0)

sigma = 0.2;
mu = 0.05;
r = 0;
dt=0.01;
S = zeros(1000,1/dt); %Initialize S
S(:,1) = 100;

rands=randn(1000,1/dt); %Create a matrix of random numbers

%W = zeros(1000,1/dt); %Initialize W 
%W(:,1) = 0;

for i=2:1/dt
    S(:,i) = S(:,i-1) + r*dt*S(:,i-1) + sigma*S(:,i-1).*rands(:,i)*sqrt(dt);
end

figure
%plot(0:dt:1-dt, S(:,1:100))
plot(0:dt:1-dt, S)

mean(S(:,end))
std(S(:, end))

%% 1. a) European Call Option

K=100;
T=1;
%logS = log(S);

% price at time 0 == T (due to r =0)
call_payoff = mean(exp(-r*T)*max(S(:, end)-K, 0))
call_BS = 100*normcdf(0.1) - 100*exp(-r*T)*normcdf(-0.1)

%% 1. b) European Asian Option

A = zeros(1000, 1);
for i = 1:size(S,1) % row
    A(i) = mean(S(i, :)); % save average price of each path
end
asian_payoff = mean(exp(-r*T)*max(A - K, 0))
