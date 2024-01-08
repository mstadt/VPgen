clear all; 

% set factor
factor = 'VII'
note = 'factorVII'

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
F_lev   = csvtable2.Var2;
F_noOC2 = csvtable3.Var2;
F_dsg   = csvtable4.Var2;

F_noOC = [F_noOC1; F_noOC2]; % merge noOC data together
allIndices = 1:length(F_noOC);

MEAN_lev = 12; %Table 1 taken from Middeldorp et al. 2000
MEAN_dsg = 32;

MEAN_err= 0.5; % difference from actual value

rng(1);

fname = '02-Jan-2024_FactorVII_newNoOC_lev_dsg.mat';

% Load best indices
best_dat = load(fname);

F_noOC_lev = best_dat.F_noOC_lev;
F_noOC_dsg = best_dat.F_noOC_dsg;

noOC_lev = sort(F_noOC_lev);
noOC_dsg = sort(F_noOC_dsg);


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

best_err = abs(mean_diff_lev - MEAN_lev) + abs(mean_diff_dsg - MEAN_dsg);
fprintf('starting best error: %0.5f \n', best_err)
% see if can improve with swapping 1 value at a time
noOC_lev_orig = noOC_lev;
noOC_dsg_orig = noOC_dsg;
for ii = 1:28
    disp(ii)

    for jj = 1:28
        % reset
        noOC_lev = noOC_lev_orig;
        noOC_dsg = noOC_dsg_orig;

        % change one VP value
        temp = noOC_lev(ii);
        noOC_lev(ii) = noOC_dsg(jj);
    
        noOC_dsg(jj) = temp;
    
        diff_lev = F_lev - noOC_lev;
        diff_dsg = F_dsg - noOC_dsg;
    
        mean_diff_lev = mean(diff_lev); % mean differences lev
        mean_diff_dsg = mean(diff_dsg); % mean differences dsg
    
        diff_lev_rd = round(mean_diff_lev, 0, 'decimals');
        diff_dsg_rd = round(mean_diff_dsg, 0, 'decimals');
    
        tot_err = abs(mean_diff_lev - MEAN_lev) + abs(mean_diff_dsg - MEAN_dsg);
    
        if tot_err < best_err
            best_err = tot_err;
            fprintf('ii: %i, jj: %i \n', ii,jj)
            fprintf('error improved, best error: %0.5f \n', best_err)
        end
    
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
    
        if and(OBJ_dsg, OBJ_lev)
            fprintf('OBJ REACHED! \n')
            break;
        end
    end

end

