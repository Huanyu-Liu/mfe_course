% shrink7.m - compare fitted slope with theoretical one

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

% regress true expected return on observed historical track record 
b2=regress(expret',[ones(n,1) averet'])
expfit2=b2(1)+b2(2)*averet;
plot(averet,expret,'.',averet,expfit2,'-',averet,averet,'--')
axis('equal')

% normalized by cross-sectional average expecter return
expret3=expret-mean(expret);
averet3=averet-mean(averet);
b3=regress(expret3',averet3')
expfit3=b3*averet3;
plot(averet3,expret3,'.',averet3,expfit3,'-',averet3,averet3,'--')
xlabel('Historical Performance Relative to Peer Average')
ylabel('True Skill Relative to Peer Average')
axis('equal')

% compare with theoretical shrinkage slope
Beta=delta^2/(delta^2+(omega^2/t))
plot(averet3,expret3,'.',averet3,expfit3,'-',averet3,averet3,'--', ...
   averet3,Beta*averet3,'-.')
xlabel('Historical Performance Relative to Peer Average')
ylabel('True Skill Relative to Peer Average')
axis('equal')

