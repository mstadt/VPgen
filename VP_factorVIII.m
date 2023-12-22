% clear
clear all
%close all

% ISSUES:
%  - Factor VIII is tough to get.
%  - Hemophilia: Factor VIII lower than 50% of normal
%  - The kernel densities definitely show hemophilia and make
%     it a high chance that this will be drawn
%  - How can we address this issue?
% -- Also --- the distributions given in Table 1 don't make sense for the
% data presented in Figure 1 (Lev before OC is 103 +/- 4 in Tab 1) but the variance
% is much larger in Figure 1
% IDEAS: remove the low patients from the data set to change the KDF....


% set factor
factor = 'VIII'
note = 'factorVIII'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data and Set Ranges %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load Data
% temp = strcat('data/Factor', factor, '_noOC1.csv');
% csvtable1 =  readtable(temp);
temp = strcat('data/Factor', factor, '_lev.csv');
csvtable2 =  readtable(temp);
% temp = strcat('data/Factor', factor, '_noOC2.csv');
% csvtable3 =  readtable(temp);
temp = strcat('data/Factor', factor, '_dsg.csv');
csvtable4 =  readtable(temp);

temp = '22-Dec-2023_FactorVIII_newNoOC_lev_dsg.mat'; % No OC sorted properly
noOCdat = load(temp);
%F_noOC1 = csvtable1.Var2; 
F_noOC1 = noOCdat.F_noOC_lev;
F_noOC2 = noOCdat.F_noOC_dsg;

F_lev   = csvtable2.Var2;
%F_noOC2 = csvtable3.Var2;
F_dsg   = csvtable4.Var2;

F_noOC = [F_noOC1; F_noOC2]; % merge noOC data together

% skew towards the mean
F_noOC = [F_noOC; mean(F_noOC); prctile(F_noOC, [60, 75, 80])'];
F_lev = [F_lev; mean(F_lev); prctile(F_lev, [60, 75])'];
F_dsg = [F_dsg; mean(F_dsg); prctile(F_dsg, [60, 75])'];

% Remove values less than 65
% F_noOC(F_noOC < 65) = [];
% F_lev(F_lev < 65) = [];
% F_dsg(F_dsg < 65) = [];


%%
%clearvars -except F_noOC F_lev F_dsg factor note;

%Desired Mean Difference After Treatment (Lev - NoOC):
% Factor II
MEAN_lev = 6; %Table 1 taken from Middeldorp et al. 2000
STD_lev  = 25;  %The Standard Deviation.

MEAN_dsg = 10;
STD_dsg  = 23;

N_vp = 1e3; %1e4; %100 %1e4; %100 % how many virtual patients

% shuffle hyperparameters
sigma1_lev = 0.1 * N_vp; %0.25 * N_vp; %
sigma2_lev = 0.05*sigma1_lev; % variance for high and low values
p_lev = [25, 75]; % percentiles to change bias

MEAN_err = 0.1; % percentage error from given mean
STD_err  = 0.1; % percentage error from given STD

sigma1_dsg = 0.1 * N_vp; %0.15 * N_vp; 
sigma2_dsg = 0.05*sigma1_dsg; 
p_dsg = [25, 75]; % percentiles to change bias



% How many trials do we want at max
MAX_TRIALS = 1e3 %1e5;

% set random seed
rng(10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create plots of original distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel density function for lev and dsg No OC
% No OC Factors 
[PDFF_noOC1,XF_noOC1] = ksdensity(F_noOC1);
[PDFF_noOC2,XF_noOC2] = ksdensity(F_noOC2);
[PDFF_noOC,XF_noOC] = ksdensity(F_noOC);

cmap = parula(6); lw = 4;
figure(50);
clf;
hold on
plot(XF_noOC1,PDFF_noOC1,'color',cmap(2,:),'linewidth',lw)
plot(XF_noOC2,PDFF_noOC2,'color',cmap(4,:),'linewidth',lw)
plot(XF_noOC, PDFF_noOC,'color',cmap(6,:),'linewidth',lw)

legend('no oc 1', 'no oc 2', 'no oc combined')


%%
% Kernel density functions
% No OC Factors 
[PDFF_noOC,XF_noOC] = ksdensity(F_noOC);

%Lev Factors
[PDFF_lev,XF_lev] = ksdensity(F_lev);

% Dsg Factors
[PDFF_dsg,XF_dsg] = ksdensity(F_dsg);


figure(1)
clf;
lw = 3;
cmap = parula(5);
c_h = 1; c_kdf = 3; c_samp = 5;
w_bin = 5;
xrange = [10, 250];
yrange = [0,0.05];
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Inverse Kernel Density Functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Range 
pi = linspace(.0001,.9999,10000);

%No OC Factor: 
xiNoOC = ksdensity(F_noOC,pi,'Function','icdf');

%Lev Factors
xiLev = ksdensity(F_lev,pi,'Function','icdf');

% Dsg Factors
xiDsg = ksdensity(F_dsg, pi, 'Function', 'icdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create random paired samples %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OBJ = 0; % Do we match objective?
NUM_TRIALS = 0; % number of trials

% Pick no OC values for each patient (ordered)
xqNoOC = rand(1, N_vp);
samplesNoOC = sort(interp1(pi,xiNoOC,xqNoOC)); % ordered no OC samples
while sum(isnan(samplesNoOC)) > 0
    %fprintf('nan value, resampling \n')
    xqNoOC = rand(1, N_vp);
    samplesNoOC = sort(interp1(pi,xiNoOC,xqNoOC));
end

%%%%%%%%%%%%%%%
% Lev samples %
%%%%%%%%%%%%%%%
% Pick lev samples 
xqLev = rand(1, N_vp);
samplesLev = sort(interp1(pi,xiLev,xqLev)); % ordered lev samples

while sum(isnan(samplesLev)) > 0
    %fprintf('nan value, resampling \n')
    xqLev = rand(1, N_vp);
    samplesLev = sort(interp1(pi,xiLev,xqLev));
end

% Compute differences
diff_lev = samplesLev - samplesNoOC;

mean_diff_lev = mean(diff_lev);
std_diff_lev  = std(diff_lev);

%%
% Plot starting mean diff and std diff
figure(2);
clf;
subplot(1,2,1)
hold on
c_p = 2;
% plot(NUM_TRIALS, mean_diff_lev, 'linestyle','none', 'marker', 'o', ...
%                             'markersize', 10,...
%                             'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
yline(MEAN_lev*(1-MEAN_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(MEAN_lev, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(MEAN_lev*(1+MEAN_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('mean difference')
title('MEAN')

subplot(1,2,2)
hold on
% plot(NUM_TRIALS, std_diff_lev, 'linestyle','none', 'marker', 'o', 'markersize', 10,...
%                                 'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
yline(STD_lev*(1-STD_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(STD_lev, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(STD_lev*(1+STD_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('std difference')
title('STD')


% Introduce variance into the lev sample
% by doing a biased shuffle for the "matching"
OBJ_MEAN = 0;
last_samp = 0;
while and(OBJ ~= 1, NUM_TRIALS < MAX_TRIALS)

    NUM_TRIALS = NUM_TRIALS + 1;
    % Prevent just shuffling if mean is good
    samp_check = NUM_TRIALS - last_samp;
    if or(OBJ_MEAN == 0, samp_check > 100)
        last_samp = NUM_TRIALS;

        % Pick no OC values for each patient (ordered)
        xqNoOC = rand(1, N_vp);
        samplesNoOC = sort(interp1(pi,xiNoOC,xqNoOC)); % ordered no OC samples
        while sum(isnan(samplesNoOC)) > 0
            %fprintf('nan value, resampling \n')
            xqNoOC = rand(1, N_vp);
            samplesNoOC = sort(interp1(pi,xiNoOC,xqNoOC));
        end
        % draw new lev and OC samples if sample did not fit mean objective
        xqLev = rand(1, N_vp);
        samplesLev = sort(interp1(pi,xiLev,xqLev)); % ordered lev samples

        while sum(isnan(samplesLev)) > 0
            %fprintf('nan value, resampling \n')
            xqLev = rand(1, N_vp);
            samplesLev = sort(interp1(pi,xiLev,xqLev));
        end
    end

    
    % shuffle the Lev samples (biased)
    samplesLev_new = biasedShuffle(samplesLev, sigma1_lev, sigma2_lev, p_lev);

    % Compute differences
    diff_lev_new = samplesLev_new - samplesNoOC; %samplesLev - samplesNoOC;
    
    mean_diff = mean(diff_lev_new);
    std_diff  = std(diff_lev_new);
    

    % Plot MEAN and STD difference
    figure(2);
    subplot(1,2,1)
    plot(NUM_TRIALS, mean_diff, 'linestyle','none', 'marker', 'o', ...
                                'markersize', 10,...
                                'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
    
    subplot(1,2,2)
    plot(NUM_TRIALS, std_diff, 'linestyle','none', 'marker', 'o', 'markersize', 10,...
                                    'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
    if mod(NUM_TRIALS,200) == 0
        pause(0.1)
    end

    % Check objectives
    OBJ_MEAN = 0; OBJ_SD = 0;
    if abs((mean_diff - MEAN_lev)/MEAN_lev) < MEAN_err 
        OBJ_MEAN = 1; 
    end
    if abs((std_diff - STD_lev)/STD_lev) < STD_err 
        OBJ_SD = 1;
    end
    
    OBJ_test = OBJ_MEAN + OBJ_SD;
    if OBJ_test == 2
        OBJ = 1;
        fprintf('objective reached in %i trials \n', NUM_TRIALS)
    end

    samplesLev = samplesLev_new; % set samplesLev to new samples
end

fprintf('lev samples complete. \n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DSG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OBJ = 0; % Do we match objective?
NUM_TRIALS = 0; % number of trials
% Pick dsg samples 
xqDsg = rand(1, N_vp);
samplesDsg = sort(interp1(pi,xiDsg,xqDsg)); % ordered dsg samples
while sum(isnan(samplesDsg)) > 0
    % resample if nan value
    %fprintf('nan value, resampling \n')
    xqDsg = rand(1, N_vp);
    samplesDsg = sort(interp1(pi,xiDsg,xqDsg));
end

% compute differences
diff_dsg = samplesDsg - samplesNoOC;

mean_diff_dsg = mean(diff_dsg);
std_dif_dsg = std(diff_dsg);

% Plot starting mean diff and std diff
figure(10);
clf;
subplot(1,2,1)
hold on
c_p = 1;
yline(MEAN_dsg*(1-MEAN_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(MEAN_dsg, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(MEAN_dsg*(1+MEAN_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('mean difference')
title('MEAN')

subplot(1,2,2)
hold on
yline(STD_dsg*(1-STD_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
yline(STD_dsg, 'color', 'black', 'linewidth', 2, 'linestyle', '-')
yline(STD_dsg*(1+STD_err), 'color', 'black', 'linewidth', 2, 'linestyle', '--')
xlabel('trial')
ylabel('std difference')
title('STD')

% Biased shuffle and redraw
OBJ_MEAN = 0;
last_samp = 0;
% TO DO: bias DSG samples based on what the lev change was for each patient
%diff_lev = samplesLev - samplesNoOC; % difference on lev for each patient
while and(OBJ ~=1, NUM_TRIALS < MAX_TRIALS)
    NUM_TRIALS = NUM_TRIALS + 1;
    % prevent just shuffling if mean is good
    samp_check = NUM_TRIALS - last_samp;
    if or(OBJ_MEAN == 0, samp_check > 200)
        last_samp = NUM_TRIALS;

        % pick new dsg values (ordered)
        xqDsg = rand(1, N_vp);
        samplesDsg = sort(interp1(pi,xiDsg,xqDsg)); % ordered dsg samples
        while sum(isnan(samplesDsg)) > 0
            % resample if nan value
            %fprintf('nan value, resampling \n')
            xqDsg = rand(1, N_vp);
            samplesDsg = sort(interp1(pi,xiDsg,xqDsg));
        end
    end

    % shuffle the dsg samples (biased)
    samplesDsg_new = biasedShuffle(samplesDsg, sigma1_dsg, sigma2_dsg, p_dsg);

    % compute differences
    diff_dsg_new = samplesDsg_new - samplesNoOC;

    mean_diff = mean(diff_dsg_new);
    std_diff = std(diff_dsg_new);

    % Plot MEAN and STD difference
    figure(10);
    subplot(1,2,1)
    plot(NUM_TRIALS, mean_diff, 'linestyle','none', 'marker', 'o', ...
                                'markersize', 10,...
                                'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
    
    subplot(1,2,2)
    plot(NUM_TRIALS, std_diff, 'linestyle','none', 'marker', 'o', 'markersize', 10,...
                                    'MarkerFaceColor', cmap(c_p,:), 'color', cmap(c_p,:))
    if mod(NUM_TRIALS,200) == 0
        pause(0.1)
    end

    % Check objectives
    OBJ_MEAN = 0; OBJ_SD = 0;
    if abs((mean_diff - MEAN_dsg)/MEAN_dsg) < MEAN_err
        OBJ_MEAN = 1;
    end
    if abs((std_diff - STD_dsg)/STD_dsg) < STD_err
        OBJ_SD = 1;
    end

    OBJ_test = OBJ_MEAN + OBJ_SD;
    if OBJ_test == 2
        OBJ = 1;
        fprintf('objective reached in %i trials \n', NUM_TRIALS)
    end

    samplesDsg = samplesDsg_new;
end
fprintf('dsg samples done. \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot samples on kernel density function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
figure(3);

diff_lev = samplesLev - samplesNoOC;
diff_dsg = samplesDsg - samplesNoOC;

yrange = [0,0.1];
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
xrange = [40, 200];
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




