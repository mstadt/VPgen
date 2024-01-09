clear all;

% set factor
factor = 'VII';
%note = 'alg3' % algorithm 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data and Set Ranges %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnoOC = './sort_noOC/02-Jan-2024_FactorVII_newNoOC_lev_dsg.mat';
dat = load(fnoOC);
NoOC_lev = dat.F_noOC_lev;
NoOC_dsg = dat.F_noOC_dsg;
NoOC = [NoOC_lev;NoOC_dsg];

[PDFF_noOC_lev,XF_noOC_lev] = ksdensity(NoOC_lev);
[PDFF_noOC_dsg,XF_noOC_dsg] = ksdensity(NoOC_dsg);
[PDFF_noOC,XF_noOC] = ksdensity(NoOC);


% plots of distributions
figure(15)
clf;
lw = 3;
cmap = parula(6);
c1 = cmap(1,:); c2 = cmap(3,:); c3 = cmap(5,:);
hold on
plot(XF_noOC_lev,PDFF_noOC_lev,'color',c1,'linewidth',lw)
plot(XF_noOC_dsg,PDFF_noOC_dsg,'color',c2,'linewidth',lw)
plot(XF_noOC,PDFF_noOC,'color',c3,'linewidth',lw)

temp = sprintf('Factor %s before OC',factor);
xlabel(temp)

legend('No OC before lev', 'No OC before dsg','Combined')
