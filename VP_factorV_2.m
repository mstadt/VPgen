% clear
clear all;

% Set factor
factor = 'V'
note = 'alg2' % algorithm 2 (based on Suzanne & Melissa meeting 2024-01-08)

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

clearvars -except F_noOC F_lev F_dsg factor note;

%Desired Mean Difference After Treatment (Lev - NoOC):
% Factor V
MEAN_lev = -3; %Table 1 taken from Middeldorp et al. 2000
STD_lev  = 12;  %The Standard Deviation.

MEAN_dsg = -11;
STD_dsg  = 8;

N_vp = 1e4; %100; %1e4; %100  % how many virtual patients

% hyperparameters
delta_lev = 0.25; % hyperparameter for added noise on percentiles
delta_dsg = 0.15; % hyperparamger

MEAN_err = 0.2; %0.1; % percentage error from given mean
STD_err  = 0.1; % percentage error from given STD

% Maximum number of trials
MAX_TRIALS = 100; %500; %50 %1e4; %100; %5e3; %100; %500; %5e3; %100;

% set random seed
rng(1) %rng(72)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute kernel density estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel density estimates
% No OC Factors 
[PDFF_noOC,XF_noOC] = ksdensity(F_noOC);

%Lev Factors
[PDFF_lev,XF_lev] = ksdensity(F_lev);

% Dsg Factors
[PDFF_dsg,XF_dsg] = ksdensity(F_dsg);

% plots of original distribution
figure(1)
clf;
lw = 3;
cmap = parula(5);
c_h = 1; c_kdf = 3; c_samp = 5;
w_bin = 5;
xrange = [40, 190];
yrange = [0,0.075];
subplot(1,3,1)
histogram(F_noOC,'Normalization','pdf', ...
            'BinWidth', w_bin, 'FaceColor', cmap(c_h,:))
hold on
plot(XF_noOC,PDFF_noOC,'color',cmap(c_kdf,:),'linewidth',lw)
xlim(xrange)
ylim(yrange)
temp = strcat('Factor ',factor, ' (No OC)');
title({'Histogram and Kernel Density Function',temp})

subplot(1,3,2)
histogram(F_lev,'Normalization','pdf',...
            'BinWidth', w_bin, 'FaceColor', cmap(c_h,:))
hold on
plot(XF_lev,PDFF_lev,'color',cmap(c_kdf,:),'linewidth',lw)
xlim(xrange)
ylim(yrange)
temp = strcat('Factor ', factor,' (Lev)');
title({'Histogram and Kernel Density Function',temp})

subplot(1,3,3)
histogram(F_dsg,'Normalization','pdf', ...
            'BinWidth', w_bin, 'FaceColor', cmap(c_h,:))
hold on
plot(XF_dsg,PDFF_dsg,'color',cmap(c_kdf,:),'linewidth',lw)
xlim(xrange)
ylim(yrange)
temp = strcat('Factor ',factor,' (Dsg)');
title({'Histogram and Kernel Density Function',temp})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create inverse kernel density functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Range 
pi = linspace(.0001,.9999,10000);

%No OC Factor: 
xiNoOC = ksdensity(F_noOC,pi,'Function','icdf');

%Lev Factors
xiLev = ksdensity(F_lev,pi,'Function','icdf');

% Dsg Factors
xiDsg = ksdensity(F_dsg, pi, 'Function', 'icdf');

%%%%%%%%%%%%%%%%%
% No OC sample
%%%%%%%%%%%%%%%%%
% Pick percentiles for No OC values
xqNoOC = lhsdesign(1,N_vp); % latin hypercube sampling %rand(1,N_vp);

% Interpolate factor levels based on xqNoOC and inverse
% kernel density function
samplesNoOC = interp1(pi, xiNoOC, xqNoOC);
while sum(isnan(samplesNoOC)) > 0
    % resample if nan value
    % fprintf('nan value, resampling \n')
    xqNoOC = lhsdesign(1,N_vp);
    samplesNoOC = interp1(pi, xiNoOC, xqNoOC);
end


%% LEV AND DSG together
%%%%%%%%%%
% LEV and DSG
%%%%%%%%%
fprintf('starting lev and dsg. \n')
% Plot mean diff and std ranges
figure(2);
clf;
subplot(1,2,1)
hold on
yline(MEAN_lev*(1-MEAN_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(MEAN_lev, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(MEAN_lev*(1+MEAN_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('mean difference (Lev - NoOC)')
title('MEAN')

subplot(1,2,2)
hold on
yline(STD_lev*(1-STD_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(STD_lev, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(STD_lev*(1+STD_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('std difference (Lev - NoOC)')
title('STD')
sgtitle('LEV')
figure(3);
clf;
subplot(1,2,1)
hold on
c_p = 2;
yline(MEAN_dsg*(1-MEAN_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(MEAN_dsg, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(MEAN_dsg*(1+MEAN_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('mean difference (Dsg - NoOC)')
title('MEAN')

subplot(1,2,2)
hold on
yline(STD_dsg*(1-STD_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(STD_dsg, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(STD_dsg*(1+STD_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('std difference (Dsg - NoOC)')
title('STD')
sgtitle('DSG')

% reset OBJ and NUMTRIALS
NUM_TRIALS = 0;
OBJ = 0;
last_samp = 0;
while OBJ~=1
    NUM_TRIALS = NUM_TRIALS + 1;
    if NUM_TRIALS > MAX_TRIALS
        fprintf('MAX TRIALS reached. \n')
        break;
    end % if MAXTRIALS

    % resample No OC if reach a certain number of trials
    samp_check = NUM_TRIALS - last_samp;
    %if samp_check > 25
    % resample no OC
    %last_samp = NUM_TRIALS;
    xqNoOC = rand(1,N_vp); %lhsdesign(1,N_vp); % latin hypercube sampling %rand(1,N_vp);
    samplesNoOC = interp1(pi, xiNoOC, xqNoOC);
    while sum(isnan(samplesNoOC)) > 0
        % resample if nan value
        % fprintf('nan value, resampling \n')
        xqNoOC = lhsdesign(1,N_vp);
        samplesNoOC = interp1(pi, xiNoOC, xqNoOC);
    end
    %end

    

    % random noise (LEV)
    % Min and max values for noise
    min_vals = max(-delta_lev, min(pi) - xqNoOC);
    max_vals = min(delta_lev, max(pi) - xqNoOC);
    DELTA_lev = min_vals + (max_vals - min_vals).* rand(1, N_vp);
    xqLev = xqNoOC + DELTA_lev; % add noise to xqNoOC

    % LEV samples
    samplesLev = interp1(pi,xiLev,xqLev); % Lev samples 
    diff_lev = samplesLev - samplesNoOC;
    mean_diff_lev = mean(diff_lev);
    std_diff_lev  = std(diff_lev);

    % random noise (DSG)
    % Min and max values for noise
    min_vals = max(-delta_dsg, min(pi) - xqNoOC);
    max_vals = min(delta_dsg, max(pi) - xqNoOC);
    DELTA_dsg = min_vals + (max_vals - min_vals).* rand(1, N_vp);
    xqDsg = xqNoOC + DELTA_dsg; % add noise to xqNoOC

    % DSG samples
    samplesDsg = interp1(pi,xiDsg,xqDsg); % dsg samples
    diff_dsg = samplesDsg - samplesNoOC;
    mean_diff_dsg = mean(diff_dsg);
    std_diff_dsg = std(diff_dsg);

    % Plot MEAN and STD diff (LEV)
    figure(2);
    subplot(1,2,1)
    plot(NUM_TRIALS, mean_diff_lev, 'linestyle','none', 'marker', 'o', ...
                                'markersize', 10,...
                                'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
    
    subplot(1,2,2)
    plot(NUM_TRIALS, std_diff_lev, 'linestyle','none', 'marker', 'o', 'markersize', 10,...
                                    'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))

    %Plot MEAN and STD difference (DSG)
    figure(3);
    subplot(1,2,1)
    plot(NUM_TRIALS, mean_diff_dsg, 'linestyle','none', 'marker', 'o', ...
                                'markersize', 10,...
                                'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
    
    subplot(1,2,2)
    plot(NUM_TRIALS, std_diff_dsg, 'linestyle','none', 'marker', 'o', 'markersize', 10,...
                                    'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
    if mod(NUM_TRIALS,300) == 0
        pause(0.1)
    end

    % Check objectives
    OBJ_MEAN_dsg = 0; OBJ_SD_dsg = 0; OBJ_dsg = 0;
    if abs((mean_diff_dsg - MEAN_dsg)/MEAN_dsg) < MEAN_err 
        OBJ_MEAN_dsg = 1; 
    end
    if abs((std_diff_dsg - STD_dsg)/STD_dsg) < STD_err 
        OBJ_SD_dsg = 1;
    end
    
    OBJ_test_dsg = OBJ_MEAN_dsg + OBJ_SD_dsg;
    if OBJ_test_dsg == 2
        OBJ_dsg = 1;
    end

    % Check objectives (lev)
    OBJ_MEAN_lev = 0; OBJ_SD_lev = 0; OBJ_lev = 0;
    if abs((mean_diff_lev - MEAN_lev)/MEAN_lev) < MEAN_err 
        OBJ_MEAN_lev = 1; 
    end
    if abs((std_diff_lev - STD_lev)/STD_lev) < STD_err 
        OBJ_SD_lev = 1;
    end
    
    OBJ_test_lev = OBJ_MEAN_lev + OBJ_SD_lev;
    if OBJ_test_lev == 2
        OBJ_lev = 1;
    end

    OBJ_test = OBJ_lev + OBJ_dsg;
    if OBJ_test == 2
        OBJ = 1;
        fprintf('objective reached in %i trials \n', NUM_TRIALS)
        fprintf('complete. \n')
    end
end % while


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot samples on kernel density function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
subplot(1,3,1)
w_bin2 = 2;
histogram(samplesNoOC,'Normalization','pdf', ...
                'BinWidth', w_bin2, 'FaceColor', cmap(c_samp,:))

subplot(1,3,2)
histogram(samplesLev, 'Normalization', 'pdf',...
                'BinWidth', w_bin2, 'FaceColor', cmap(c_samp,:))

subplot(1,3,3)
histogram(samplesDsg, 'Normalization', 'pdf', ...
                'BinWidth', w_bin2, 'FaceColor', cmap(c_samp,:))

legend({'Data', 'Kernel Density function', 'VP'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot differences for pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure(7);

diff_lev = samplesLev - samplesNoOC;
diff_dsg = samplesDsg - samplesNoOC;

yrange = [0,0.075];
clf; 
subplot(1,2,1)
histogram(diff_lev, ...
                'BinWidth', w_bin2, 'FaceColor', cmap(2,:), ...
                'Normalization', 'pdf')

xlabel(strcat('Factor ', factor,' level difference after lev'))
ylabel('frequency')
ylim(yrange)
temp = sprintf('LEV \n MEAN DIFF: %0.3f \n STD DIFF: %0.3f',...
        mean(diff_lev), std(diff_lev));
title(temp)

subplot(1,2,2)
histogram(diff_dsg, ...
                'BinWidth', w_bin2, 'FaceColor', cmap(2,:),...
                'Normalization', 'pdf')
xlabel('factor level difference after dsg')
ylabel('density')
ylim(yrange)
temp = sprintf('DSG \n MEAN DIFF: %0.3f \n STD DIFF: %0.3f',...
        mean(diff_dsg), std(diff_dsg));
title(temp)

figure(13)
clf; 
hold on
histogram(diff_lev, ...
                'BinWidth', w_bin2, 'FaceColor', cmap(2,:), ...
                'Normalization', 'pdf')
histogram(diff_dsg, ...
                'BinWidth', w_bin2, 'FaceColor', cmap(4,:),...
                'Normalization', 'pdf')

xlabel(strcat('Factor ', factor,' level difference after OC'))
ylabel('density')
legend('Lev', 'Dsg')
ylim(yrange)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure(5)
clf;
subplot(1,2,1)
xrange = [30, 200];
yrange = xrange;
plot(samplesNoOC, samplesLev, 'linestyle', 'none', 'marker', '.', 'markersize', 15)
xlabel(strcat('Factor ', factor,' before OC'))
ylabel(strcat('Factor ',factor,' after Lev'))
title('Lev VP pairs')
ylim(yrange)
xlim(xrange)

subplot(1,2,2)
plot(samplesNoOC, samplesDsg, 'linestyle', 'none', 'marker', '.', 'markersize', 15)
xlabel(strcat('Factor ',factor,' before OC'))
ylabel(strcat('Factor ',factor,' after Dsg'))
title('Dsg VP pairs')
ylim(yrange)
xlim(xrange)


figure(6)
clf;
hold on
plot(samplesNoOC, samplesLev, 'linestyle', 'none', ...
    'marker', '.', 'markersize', 15, ...
    'color', cmap(2,:))
plot(samplesNoOC, samplesDsg, 'linestyle', 'none', ...
    'marker', '.', 'markersize', 15, ...
    'color', cmap(4,:))
xlabel(strcat('Factor ', factor,' before OC'))
ylabel(strcat('Factor ',factor,' after OC'))
title({'VP pairs', ['Factor ', factor]})
legend('Lev','Dsg') 
ylim(yrange)
xlim(xrange)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
save_file = 1;
fname = strcat(date, '_Factor', factor, '_VP', '_n-', ...
                num2str(N_vp), '_note-', note, '.mat');
if isfile(fname)
    save_file = input('file exists. save file? (0/1)');
    if save_file
        note = input('change note. note: ');
        fname = strcat(date, '_Factor', factor, '_VP', '_n-', ...
                    num2str(N_vp), '_note-', note, '.mat');
    end
end
if save_file
    save(fname, 'samplesNoOC', 'samplesDsg', 'samplesLev')
    fnameinf = strcat(date, '_Factor', factor, '_VP', ...
        '_n-', num2str(N_vp), '_note-', note ,'_info.mat');
    save(fnameinf)
    
    fprintf('VP saved to: \n %s \n', fname)
end
