% can I even get a NoOC combination with a Lev combination that
% matches the mean difference reported in Middeldorp et al?
clear all; 

% set factor
factor = 'VIII'
note = 'factorVIII'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data and Set Ranges %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load Data
temp = strcat('data/Factor', factor, '_noOC1.csv');
csvtable1 =  readtable(temp);
temp = strcat('data/Factor', factor, '_lev.csv');
csvtable2 =  readtable(temp);
temp = strcat('data/Factor', factor, '_noOC2.csv');
csvtable3 =  readtable(temp);
temp = strcat('data/Factor', factor, '_dsg.csv');
csvtable4 =  readtable(temp);
F_noOC1 = csvtable1.Var2; 
F_lev   = csvtable2.Var2;
F_noOC2 = csvtable3.Var2;
F_dsg   = csvtable4.Var2;

F_noOC = [F_noOC1; F_noOC2]; % merge noOC data together


N = length(F_lev);

MEAN_lev = 6; %Table 1 taken from Middeldorp et al. 2000
STD_lev  = 25;  %The Standard Deviation.

MEAN_dsg = 10;
STD_dsg  = 23;

MEAN_err= 0.05;

rng(1);

MAX_TRIALS = 1e3;

OBJ = 0;

NUM_TRIALS = 0;

cmap = parula(10);
c_p = 3;
figure(2);
clf;
hold on
% plot(NUM_TRIALS, mean_diff_lev, 'linestyle','none', 'marker', 'o', ...
%                             'markersize', 10,...
%                             'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
yline(MEAN_lev*(1-MEAN_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(MEAN_lev, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(MEAN_lev*(1+MEAN_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('mean difference')
title('MEAN')

while and(NUM_TRIALS < MAX_TRIALS, OBJ ~= 1)
    NUM_TRIALS = NUM_TRIALS + 1;

    randIndices = randperm(length(F_noOC), N); 
    
    noOC = F_noOC(randIndices);

    diff_lev = F_lev - noOC;

    mean_diff = mean(diff_lev);

    % Plot MEAN and STD difference
    figure(2);
    plot(NUM_TRIALS, mean_diff, 'linestyle','none', 'marker', 'o', ...
                                'markersize', 10,...
                                'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
    

    % check objective
    if abs((mean_diff - MEAN_lev)/MEAN_lev) < MEAN_err
        OBJ = 1;
        fprintf('objective reached in %i trials \n', NUM_TRIALS)
    end
end

F_noOC_lev = noOC;

allIndices = 1:length(F_noOC);
notSelectedIndices = setdiff(allIndices, randIndices);

F_noOC_dsg = F_noOC(notSelectedIndices);

OBJ = 0;
NUM_TRIALS = 0;
figure(3);
clf;
hold on
c_p = 7;
% plot(NUM_TRIALS, mean_diff_lev, 'linestyle','none', 'marker', 'o', ...
%                             'markersize', 10,...
%                             'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
yline(MEAN_dsg*(1-MEAN_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(MEAN_dsg, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(MEAN_dsg*(1+MEAN_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('mean difference')
title('MEAN')
%while and(NUM_TRIALS < MAX_TRIALS, OBJ ~= 1)
for Ind = 1:length(F_noOC_dsg)
    NUM_TRIALS = NUM_TRIALS + 1;

    %randInd = randi(length(F_noOC_dsg));

    noOC_dsg = F_noOC_dsg;  % copy of F_noOC_dsg
    noOC_dsg(Ind) = [];

    diff_dsg = F_dsg - noOC_dsg;
    mean_diff = mean(diff_dsg);

    figure(3);
    plot(NUM_TRIALS, mean_diff, 'linestyle','none', 'marker', 'o', ...
                                'markersize', 10,...
                                'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))

    % check objective
    if abs((mean_diff - MEAN_dsg)/MEAN_dsg) < MEAN_err
        OBJ = 1;
        fprintf('objective reached in %i trials \n', NUM_TRIALS)
        break;
    end

end

F_noOC_dsg = noOC_dsg;
%F_noOC_lev = F_noOC_lev;

fname = strcat(date, '_FactorVIII_newNoOC_lev_dsg.mat');
save(fname, 'F_noOC_lev', 'F_noOC_dsg')
fprintf('%s \n',fname)
