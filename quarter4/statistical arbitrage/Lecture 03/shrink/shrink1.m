% shrink1.m - generate cross-section of portfolio manager talent

% parameters
n=500;      % number of portfolio managers
t=60;       % number of months of track record (5 years)
mu=0.5;     % percentage return per month for average manager (6% per year)
delta=0.25; % cross-sectional standard deviation of skill levels
omega=2;    % standard deviation of monthly returns (~7% per year)

% generate cross-section of expected returns
expret=mu+delta*randn(1,n);
histogram(expret)
xlabel('Expected Return = True Manager Skill (percent per month)')
disp([mean(expret) std(expret)])

