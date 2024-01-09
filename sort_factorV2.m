% Try to get pairings of pre-lev and pre-dsg that matches stats reported
% in Tab 1

clear all; 

% set factor
factor = 'V'
note = 'factorV'

rng(1);

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
F_lev   = sort(csvtable2.Var2);
F_noOC2 = csvtable3.Var2;
F_dsg   = sort(csvtable4.Var2);

F_noOC = [F_noOC1; F_noOC2]; % merge noOC data together


MEAN_beforeOC_lev = 105; % Table 1 Lev before OC
STD_beforeOC_lev = 21; % Table 1 Lev before OC

MAX_TRIALS = 1e3;

max_err = 0.1;

figure(1);
clf;
hold on;
xlabel('trial')
ylabel('error (mean)')
title('mean error')

figure(2);
clf;
hold on;
xlabel('trial')
ylabel('error (std)')
title('std error')

cmap = parula(10);
cp = 7;

for NUM_TRIALS = 1:MAX_TRIALS
    
    randIndices = randperm(length(F_noOC), length(F_lev));
    noOC_lev = F_noOC(randIndices);

    beforeOC_mean = mean(noOC_lev);
    beforeOC_std = std(noOC_lev);

    err_mean =beforeOC_mean - MEAN_beforeOC_lev;
    err_std =beforeOC_std - STD_beforeOC_lev;

    figure(1);
    plot(NUM_TRIALS, err_mean, 'linestyle', 'none', 'marker','o',...
                    'markersize',10,...
                    'markerfacecolor',cmap(cp,:), 'color',cmap(cp,:))

    figure(2);
    plot(NUM_TRIALS, err_std, 'linestyle', 'none', 'marker','o',...
                    'markersize',10,...
                    'markerfacecolor',cmap(cp,:), 'color',cmap(cp,:))
    if abs(err_mean) < max_err
        %fprintf('obj reached! \n')
        %Inds_final = randIndices;
        OBJ_mean = 1;
        break;
    else
        OBJ_mean = 0;
    end

%     if abs(err_std) < max_err
%         OBJ_std = 1;
%     else
%         OBJ_std = 0;
%     end

%     if and(OBJ_mean, OBJ_std)
%         fprintf('OBJ reached! \n')
%         break;
%     end

    
end
allIndices = 1:length(F_noOC);
notSelectedIndices = setdiff(allIndices, randIndices);
noOC_dsg = sort(F_noOC(notSelectedIndices));