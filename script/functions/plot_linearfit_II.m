function [h_l h_t m b sm sb r2] = plot_linearfit_II(varargin)
% [h_l h_t m b sm sb r2] = plot_linearfit(X,y,pm,pb,pr)
% Calls function lsqfitma.m to find type II linear fit line and plots the 
% line along with its equation (with precisions pm and pb) and r^2 value 
% (precision pr). Data are normalized by a measure of their variance before
% type II linear regression.
% 
% h_l = line handle
% h_t = vector containing the three text handles
% m = slope of regression
% b = intercept of regression
% sm = 2 standard deviations of the slope
% sb = 2 standard deviations of the intercept
% r2 = r^2 of regression

g = ~isnan(varargin{1}) & ~isnan(varargin{2});
X = varargin{1}(g);
y = varargin{2}(g);

switch nargin
    case 2
        pm = 3;
        pb = 3;
        pr = 2;
    case 3
        pm = varargin{3};
        pb = varargin{3};
        pr = 2;
    case 4
        pm = varargin{3};
        pb = varargin(4);
        pr = 2;
    case 5
        pm = varargin{3};
        pb = varargin{4};
        pr = varargin{5};
end
    

if size(X,1) < size(X,2)
    X = X';
end
if size(y,1) < size(y,2)
    y = y';
end

g = ~isnan(X) & ~isnan(y) & ~isinf(X) & ~isinf(y);
X = X(g); y = y(g);

xvar = std(X);
yvar = std(y);
[m,b,r,sm,sb,~,~] = lsqfitma(X/xvar,y/yvar);
sm = sm * 2; sb = sb * 2;
m = m*yvar/xvar;
sm = sm*yvar/xvar;
b = b*yvar;
sb = sb*yvar;

r2 = r^2;
% [b, bint, ~, ~, stats] = regress(y,[ones(size(y)) X]);
x_new = [min(X) max(X)];
y_new = b + m*x_new;
h_l = plot(x_new,y_new,'k');
brange = sb;
mrange = sm;

%% ELABORATE CODE TO FORCE EVERYTHING TO THE RIGHT DECIMAL
mlog = floor(log10(abs(m)));
mshift =  pm-mlog-1;
mround = 10^-mshift*(round(m*10^mshift));
newp = mshift*(mshift>0);
mstr = num2str(mround,['%.' num2str(newp) 'f']);

mrangeround = 10^-mshift*(round(mrange*10^mshift));
mrangestr = num2str(mrangeround,['%.' num2str(newp) 'f']);

blog = floor(log10(abs(b)));
bshift =  pb-blog-1;
bround = 10^-bshift*(round(b*10^bshift));
newp = bshift*(bshift>0);
if bround < 0
    plusminus = ' - ';
    bround = -bround;
else
    plusminus = ' + ';
end
bstr = num2str(bround,['%.' num2str(newp) 'f']);

brangeround = 10^-bshift*(round(brange*10^bshift));
brangestr = num2str(brangeround,['%.' num2str(newp) 'f']);
rstr = num2str(r2,['%.' num2str(pr) 'f']);

% PRINT EQUATION ON FIGURE
h_t(1) = text(0.05,0.95,['y = ' mstr '±' mrangestr ...
    'x' plusminus bstr '±' brangestr],'sc');
% h_t(1) = text(0.05,0.95,['y = ' mstr ...
%     'x' plusminus bstr],'sc');
% h_t(1) = text(0.05,0.95,['m = ' mstr '±' mrangestr],'sc');
% PRINT R^2 ON FIGURE
h_t(2) = text(0.05,0.85,['r^2 = ' rstr],'sc');

% % PRINT P VALUE ON FIGURE
% h_t(3) = text(0.55,0.85,['p = ' num2str(stats(3),2)],'sc');

set(h_t,'fontsize',11)
set(h_l,'linewidth',3)

