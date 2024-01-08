% Get no OC combination that matches the means reported in Middeldorp et
% al. 2000 Table 1 for Factor VII
% NOTE: means will be rounded to nearest whole number as done in Table 1
% for both LEV and DSG
clear all; 

% set factor
factor = 'X'
note = 'factorX'

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
allIndices = 1:length(F_noOC);

MEAN_lev = 22; %Table 1 taken from Middeldorp et al. 2000
MEAN_dsg = 25;

MEAN_err= 0.5; % difference from actual value

rng(1);

MAX_TRIALS =1e3; %1e3; % 3e3;


NUM_TRIALS = 0;

cmap = parula(10);
% % MEAN figure
% figure(2);
% clf;
% hold on
% yline(MEAN_lev-MEAN_err, 'color', 'black', 'linewidth', 2, 'linestyle', '--')
% yline(MEAN_lev, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
% yline(MEAN_lev+MEAN_err, 'color', 'black', 'linewidth', 2, 'linestyle', '--')
% xlabel('trial')
% ylabel('mean difference')
% title('MEAN (lev)')
% 
% % DSG figure
% figure(3);
% clf;
% hold on
% yline(MEAN_dsg-MEAN_err, 'color', 'black', 'linewidth', 2, 'linestyle', '--')
% yline(MEAN_dsg, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
% yline(MEAN_dsg+MEAN_err, 'color', 'black', 'linewidth', 2, 'linestyle', '--')
% xlabel('trial')
% ylabel('mean difference')
% title('MEAN (dsg)')


% Err figure
figure(4);
clf;
hold on
yline(-MEAN_err, 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(0, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(MEAN_err, 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('mean difference')
title('Error ')


c_p1 = 7;
c_p2 = 3;
ms = 'o';

OBJ = 0;
OBJ_dsg = 0;
OBJ_lev = 0;
best_err = 100;
while OBJ~= 1
    NUM_TRIALS = NUM_TRIALS + 1;
    if NUM_TRIALS > MAX_TRIALS
        fprintf('max trials reached! \n')
        fprintf('best error: %0.3f \n', best_err)

        % do best noOC break
        noOC_lev = best_lev;
        noOC_dsg = best_dsg;

        break;
    end

    randIndices = randperm(length(F_noOC), length(F_lev));
    notSelectedIndices = setdiff(allIndices, randIndices);

    noOC_lev = sort(F_noOC(randIndices));
    noOC_dsg = sort(F_noOC(notSelectedIndices));

    diff_lev = F_lev - noOC_lev;
    diff_dsg = F_dsg - noOC_dsg;

    mean_diff_lev = mean(diff_lev); % mean differences lev
    mean_diff_dsg = mean(diff_dsg); % mean differences dsg


    diff_lev_rd = round(mean_diff_lev, 0, 'decimals');
    diff_dsg_rd = round(mean_diff_dsg, 0, 'decimals');

    % Check objective
    if abs(diff_lev_rd - MEAN_lev) < 0.3
        OBJ_lev = 1;
    else
        OBJ_lev = 0;
    end

    if abs(diff_dsg_rd - MEAN_dsg) < 0.3
        OBJ_dsg = 1;
    else 
        OBJ_dsg = 0;
    end


    % check if combo is better than best
    tot_err = abs(mean_diff_lev - MEAN_lev) + abs(mean_diff_dsg - MEAN_dsg);
    if tot_err < best_err
        best_inds = randIndices;
        best_lev = noOC_lev;
        best_dsg = noOC_dsg;
        best_err = tot_err;
        ms = '^';
    end

    figure(4);
    plot(NUM_TRIALS, mean_diff_dsg-MEAN_dsg, 'linestyle', 'none', 'marker', ms,...
                                'markersize', 10,...
                                'markerfacecolor', cmap(c_p1,:),'color',cmap(c_p1,:))
    plot(NUM_TRIALS, mean_diff_lev-MEAN_lev, 'linestyle', 'none', 'marker', ms,...
                                'markersize', 10,...
                                'markerfacecolor', cmap(c_p2,:),'color',cmap(c_p2,:))
    ms = 'o';



    if and(OBJ_lev, OBJ_dsg)
        OBJ = 1;
    else
        OBJ = 0;
    end

    

end

F_noOC_lev = noOC_lev;

F_noOC_dsg = noOC_dsg;

mean_diff_lev = mean(F_lev - F_noOC_lev);
mean_diff_dsg = mean(F_dsg - F_noOC_dsg);
fprintf('mean diff lev: %0.2f \n', mean_diff_lev)
fprintf('mean diff dsg: %0.2f \n', mean_diff_dsg)

tot_err = abs(mean_diff_lev - MEAN_lev) + abs(mean_diff_dsg - MEAN_dsg);
fprintf('total error: %0.2f \n', tot_err)


fname = strcat(date, '_Factor',factor,'_newNoOC_lev_dsg.mat');
save(fname, 'F_noOC_lev', 'F_noOC_dsg')
fprintf('%s \n',fname)
