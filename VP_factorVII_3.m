% clear
clear all;

% set factor
factor = 'VII'
note = 'alg3' % algorithm 3

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
% Factor VII
MEAN_lev = 12; %Table 1 taken from Middeldorp et al. 2000
MEAN_lev_range = [MEAN_lev - 1.3, MEAN_lev +  1.0];
STD_lev  = 15;  %The Standard Deviation.
STD_lev_range = [STD_lev - 1.0, STD_lev + 1.0;];

MEAN_dsg = 32;
MEAN_dsg_range = [MEAN_dsg - 1.0, MEAN_dsg + 1.0];
STD_dsg  = 10;
STD_dsg_range = [STD_dsg - 1.0, STD_dsg + 1.0;];

N_vp = 1e4; %100; %1e4; %100; %1e4; %100  % how many virtual patients

% hyperparameters
delta_lev = 0.1; % hyperparameter for added noise on probabilities
delta_dsg = 0.1; % hyperparamger


% Maximum number of trials
MAX_TRIALS = 100; %500;

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
temp = sprintf('Factor %s before OC',factor);
xlabel(temp)
%title({'Histogram and Kernel Density Function',temp})

subplot(1,3,2)
histogram(F_lev,'Normalization','pdf',...
            'BinWidth', w_bin, 'FaceColor', cmap(c_h,:))
hold on
plot(XF_lev,PDFF_lev,'color',cmap(c_kdf,:),'linewidth',lw)
xlim(xrange)
ylim(yrange)
temp = sprintf('Factor %s after Lev',factor);
xlabel(temp)
%title({'Histogram and Kernel Density Function',temp})

subplot(1,3,3)
histogram(F_dsg,'Normalization','pdf', ...
            'BinWidth', w_bin, 'FaceColor', cmap(c_h,:))
hold on
plot(XF_dsg,PDFF_dsg,'color',cmap(c_kdf,:),'linewidth',lw)
xlim(xrange)
ylim(yrange)
temp = sprintf('Factor %s after Dsg',factor);
xlabel(temp)

sgtitle({'Histogram and Kernel Density Function',sprintf('Factor %s', factor)})

%%
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

% Plot inverse kernel density
figure(2)
lw=2;
clf;
subplot(1,3,1)
hold on
plot(pi, xiNoOC, 'color', 'black', 'linewidth',lw)
temp = sprintf('Factor %s (No OC)',factor);
ylabel(temp)

subplot(1,3,2)
hold on
plot(pi,xiLev, 'color', 'black', 'linewidth',lw)
temp = sprintf('Factor %s (Lev)',factor);
ylabel(temp)

subplot(1,3,3)
hold on
plot(pi,xiDsg, 'color', 'black', 'linewidth',lw)
temp = sprintf('Factor %s (Dsg)',factor);
ylabel(temp)
sgtitle('Inverse CDF')

%%
%%%%%%%%%%%%%%%%%%%
% Get no OC sample
%%%%%%%%%%%%%%%%%%%
fprintf('starting MEAN objective. \n')
% Plot MEAN ranges
figure(3);
c_p = 2;
clf;
subplot(1,2,1)
hold on
yline(MEAN_lev_range(1), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(MEAN_lev, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(MEAN_lev_range(2), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('mean difference (Lev - NoOC)')
title('MEAN')

subplot(1,2,2)
hold on
yline(MEAN_dsg_range(1), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(MEAN_dsg, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(MEAN_dsg_range(2), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('mean difference (Dsg - NoOC)')
title('MEAN')

% Objective: find xqNoOC that satisfy mean difference
OBJ_mean = 0;
NUM_TRIALS = 0;
flag = 0;
fprintf('start OBJ_mean.\n')
while ~OBJ_mean
    NUM_TRIALS = NUM_TRIALS + 1;
    if NUM_TRIALS > MAX_TRIALS
        fprintf('MAX TRIALS reached. \n')
        break;
    end % if MAXTRIALS

    % Pick percentiles for No OC values
    xqNoOC = rand(1,N_vp);

    % Interpolate factor levels based on xqNoOC and inverse
    % kernel density function
    samplesNoOC = interp1(pi, xiNoOC, xqNoOC);
    while sum(isnan(samplesNoOC)) > 0
        % resample nan values
        % fprintf('nan value, resampling \n')
        nan_ids = find(isnan(samplesNoOC));
        if length(nan_ids) > 10
            error('nanids longer than 10.\n')
        end


        for ii = 1:length(nan_ids)
            id = nan_ids(ii);
            xqNoOC(id) = rand(1); % new value
        end
        samplesNoOC = interp1(pi, xiNoOC, xqNoOC);
    end

    samplesLev = interp1(pi,xiLev,xqNoOC);
    mean_diff_lev = mean(samplesLev - samplesNoOC);

    samplesDsg = interp1(pi,xiDsg,xqNoOC);
    mean_diff_dsg = mean(samplesDsg - samplesNoOC);

    % Plot MEANS
    figure(3);
    subplot(1,2,1)
    plot(NUM_TRIALS, mean_diff_lev, 'linestyle', 'none', 'marker', 'o',...
                            'markersize', 10,...
                            'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
    subplot(1,2,2)
    plot(NUM_TRIALS, mean_diff_dsg, 'linestyle', 'none', 'marker', 'o',...
                            'markersize', 10,...
                            'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
    
    OBJ_lev = 0; OBJ_dsg = 0;
    if and(MEAN_lev_range(1) < mean_diff_lev , mean_diff_lev < MEAN_lev_range(2))
        OBJ_lev = 1;
    end
    if and(MEAN_dsg_range(1) < mean_diff_dsg , mean_diff_dsg < MEAN_dsg_range(2))
        OBJ_dsg = 1;
    end

    if and(OBJ_lev, OBJ_dsg)
        OBJ_mean = 1;
        fprintf('objective reached in %i trials \n', NUM_TRIALS)
        fprintf('OBJ_mean complete \n')
    end
end % while ~OBJ_mean

%%%%%%%%%%%%%%%%%%%%
% Pair matching
%%%%%%%%%%%%%%%%%%%%
% Objective add noise to xqNoOC so that standard deviation
%   in the differences matches Middeldorp et al., 2000 report

% Plot STD ranges
figure(4);
c_p = 2;
clf;
subplot(1,2,1)
hold on
yline(STD_lev_range(1), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(STD_lev, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(STD_lev_range(2), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('std difference (Lev - NoOC)')
title('STD')

subplot(1,2,2)
hold on
yline(STD_dsg_range(1), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(STD_dsg, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(STD_dsg_range(2), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('STD difference (Dsg - NoOC)')
title('STD')

NUM_TRIALS = 0;
OBJ_std = 0;
fprintf('start OBJstd \n')
while OBJ_std~=1
    NUM_TRIALS = NUM_TRIALS + 1;
    if NUM_TRIALS > MAX_TRIALS
        fprintf('MAX TRIALS reached. \n')
        break;
    end % if MAXTRIALS


    % random noise (LEV)
    % min and max values for noise
    min_vals = max(-delta_lev, min(pi) - xqNoOC);
    max_vals = min(delta_lev, max(pi) - xqNoOC);
    DELTA_lev = min_vals + (max_vals - min_vals).* rand(1, N_vp);
    xqLev = xqNoOC + DELTA_lev; % add noise to xqNoOC

    samplesLev = interp1(pi,xiLev,xqLev); % Lev samples
    diff_lev = samplesLev - samplesNoOC;
    std_diff_lev = std(diff_lev);

    % random noise (DSG)
    % min and max values for noise
    min_vals = max(-delta_dsg, min(pi) - xqNoOC);
    max_vals = min(delta_dsg, max(pi) - xqNoOC);
    DELTA_dsg = min_vals + (max_vals - min_vals).* rand(1, N_vp);
    xqDsg = xqNoOC + DELTA_dsg; % add noise to xqNoOC

    % DSG samples
    samplesDsg = interp1(pi,xiDsg,xqDsg); % dsg samples
    diff_dsg = samplesDsg - samplesNoOC;
    std_diff_dsg = std(diff_dsg);

    % Plot STD diff
    figure(4);
    subplot(1,2,1)
    plot(NUM_TRIALS, std_diff_lev, 'linestyle', 'none', 'marker', 'o',...
                    'markersize', 10,...
                    'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))

    subplot(1,2,2)
    plot(NUM_TRIALS, std_diff_dsg,'linestyle', 'none', 'marker','o',...
                    'markersize', 10,...
                    'markerfacecolor', cmap(c_p,:), 'color', cmap(c_p,:))

    % Check objectives
    OBJ_lev = 0; OBJ_dsg = 0;
    if and(STD_lev_range(1) < std_diff_lev, std_diff_lev < STD_lev_range(2))
        OBJ_lev = 1;
    end
    if and(STD_dsg_range(1) < std_diff_dsg, std_diff_dsg < STD_dsg_range(2))
        OBJ_dsg = 1;
    end

    if and(OBJ_lev, OBJ_dsg)
        OBJ_std = 1;
        fprintf('objective reached in %i trials \n', NUM_TRIALS)
        fprintf('OBJ_std complete \n')
    end
end % while OBJ_std ~=1


%% Plot results
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

%%%%%%%%%%%%%%%
% Plot on iCDF
%%%%%%%%%%%%%%%
cp = 2;
figure(2);
subplot(1,3,1)
plot(xqNoOC, samplesNoOC, 'linestyle', 'none','marker','o',...
                            'markersize', 5,...
                            'MarkerFacecolor', cmap(cp,:),...
                            'color', cmap(cp,:))
subplot(1,3,2)
plot(xqLev, samplesLev, 'linestyle', 'none','marker','o',...
                            'MarkerFacecolor', cmap(cp,:),...
                            'color', cmap(cp,:))
subplot(1,3,3)
plot(xqDsg, samplesDsg, 'linestyle', 'none','marker','o',...
                            'MarkerFacecolor', cmap(cp,:),...
                            'color', cmap(cp,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot differences for pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure(7);

diff_lev = samplesLev - samplesNoOC;
diff_dsg = samplesDsg - samplesNoOC;

yrange = [0,0.04];
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

xlabel(sprintf('Factor %s level difference after OC',factor))
ylabel('density')
legend('Lev', 'Dsg')
ylim(yrange)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xrange = [40, 200];
yrange = xrange;
figure(6)
clf;
hold on
x = linspace(xrange(1),xrange(2), 100);
y = x;
temp = gray(3);
cgray = temp(2,:);
plot(x,y, 'linewidth',1.0, 'color',cgray)
plot(samplesNoOC, samplesLev, 'linestyle', 'none', ...
    'marker', '.', 'markersize', 15, ...
    'color', cmap(2,:))
plot(samplesNoOC, samplesDsg, 'linestyle', 'none', ...
    'marker', '.', 'markersize', 15, ...
    'color', cmap(4,:))
xlabel(sprintf('Factor %s before OC',factor))
ylabel(sprintf('Factor %s after OC',factor))
title({'VP pairs', ['Factor ', factor]})
legend('','Lev','Dsg') 
ylim(yrange)
xlim(xrange)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if and(OBJ_mean, OBJ_std)
    save_file = 1;
else
    save_file = 0;
end
if save_file
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
end


