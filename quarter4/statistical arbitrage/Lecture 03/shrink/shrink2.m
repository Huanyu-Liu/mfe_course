% shrink2.m - generate observed monthly performances

% parameters
n=500;      % number of portfolio managers
t=60;       % number of months of track record (3 years)
mu=0.5;     % percentage return per month for average manager (6% per year)
delta=0.25; % cross-sectional standard deviation of skill levels
omega=2;    % standard deviation of monthly returns (~7% per year)

% generate cross-section of expected returns
expret=mu+delta*randn(1,n);
histogram(expret)
xlabel('Expected Return = True Manager Skill (percent per month)')
disp([mean(expret) std(expret)])

% generate realized monthly returns
realret=repmat(expret,[t 1])+omega*randn(t,n);
histogram(realret(:))
xlabel('Observed Monthly Returns (percent per month)')
disp([mean(realret(:)) std(realret(:))])
