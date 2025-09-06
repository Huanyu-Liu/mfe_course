% random student generator
%
% for UCLA Anderson MFE

clear
clf

% initialize
dirname='Roster';
rng('default')
rng('shuffle')
d=dir([dirname '\*.jpg']);
nd=length(d);
figure(1)
axis off
axis image
pause

% infinite loop
while 1
    ind=randperm(nd);
    for i=ind
        a=imread([dirname '\' d(i).name]);
        imagesc(a)
        axis off
        axis image
        t=d(i).name(1:end-9);
        t(t=='_')=' ';
        th=title(t);
        set(th,'fontsize',24)
        pause
    end
end
    
    
