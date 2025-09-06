% shrink4.m - regress observed average return on true expected return

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

% compute average realized monthly returns
averet=mean(realret);
histogram(averet)
xlabel('Average Realized Return = Track Record (percent per month)')
disp([mean(averet)  std(averet)])

% regress average realized return on expected return
b1=regress(averet',[ones(n,1) expret'])
avefit1=b1(1)+b1(2)*expret;
plot(expret,averet,'.',expret,avefit1,'-')
xlabel('True Skill')
ylabel('Historical Performance')
axis('equal')
