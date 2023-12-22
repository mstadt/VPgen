% Generate Factor VP samples
clear all
close all
%test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  User Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FactorName = input('what factor level? (II, V, VII, VIII, X):') %'VIII'; % II, V, VII, VIII, X
N_VP       = 100; % how many patients

% tolerance
TOL_lev     = 0.1; %0.1;
TOL_dsg     = 0.1; %0.1;
do_plots = 1; % 0/1 if want plots

% Max trials we want
MAX_TRIALS    = 100 %1e6; %1e5; %4e6; %5e6;

per_mutants_lev = 0.50; % percentatge of VPs mutated 
per_mutants_dsg = 0.70; % percentages of VPs mutated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mut_cnt = 0; % mutation count
save_VP = 1;

% get the Middeldorp et al 2000 Mean and SD difference after treatment
% (lev/dsg - noOC)
if strcmp(FactorName,'II')
    MEAN_DIFF_LEV = 12; % Table 1, FII, lev
    SD_DIFF_LEV   = 8;  % Table 1, FII, lev

    MEAN_DIFF_DSG = 16; % Table 1, FII, dsg
    SD_DIFF_DSG   = 6;  % Table 1, FII, dsg
elseif strcmp(FactorName, 'V')
    MEAN_DIFF_LEV = -3; % Table 1, FV, lev
    SD_DIFF_LEV   = 12; % Table 1, FV, lev

    MEAN_DIFF_DSG = -11; % Table 1, FV, dsg
    SD_DIFF_DSG   = 8;   % Table 1, FV, dsg
elseif strcmp(FactorName, 'VII')
    MEAN_DIFF_LEV = 12; % Table 1, FVII, lev
    SD_DIFF_LEV   = 15; % Table 1, FVII, lev

    MEAN_DIFF_DSG = 32; % Table 1, FVII, dsg
    SD_DIFF_DSG   = 10; % Table 1, FVII, dsg
elseif strcmp(FactorName,'VIII')
    MEAN_DIFF_LEV = 6; % Table 1, FVIII, lev
    SD_DIFF_LEV   = 25; % Table 1, FVIII, lev
    
    MEAN_DIFF_DSG = 10; % Table 1, FVIII, dsg
    SD_DIFF_DSG   = 23;   % Table 1, FVIII, dsg
elseif strcmp(FactorName, 'X')
    MEAN_DIFF_LEV = 22; % Table 1, FX, lev
    SD_DIFF_LEV   = 14; % Table 1, FX, lev

    MEAN_DIFF_DSG = 25; % Table 1, FX, dsg
    SD_DIFF_DSG   = 12; % Table 1, FX, dsg
else
    error('FactorName not specified')
end


% Set the random seed
rng(16);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data and Set Ranges
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname1 = strcat('./data/Factor',FactorName,'_noOC1.csv');
csvtable1 = readtable(fname1);
fname2 = strcat('./data/Factor',FactorName,'_lev.csv');
csvtable2 = readtable(fname2);
fname3 = strcat('./data/Factor',FactorName,'_noOC2.csv');
csvtable3 = readtable(fname3);
fname4 = strcat('./data/Factor',FactorName,'_dsg.csv');
csvtable4 = readtable(fname4);
F_noOC1 = csvtable1.Var2;
F_lev   = csvtable2.Var2;
F_noOC2 = csvtable3.Var2;
F_dsg   = csvtable4.Var2;

% Put no OC samples into one vector
F_noOC = [F_noOC1; F_noOC2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create plots of original distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot Kernel Density Function
if do_plots
    fprintf('** plotting kernel density functions and data ** \n')
    xrange = [40, 190];
    yrange = [0, 0.03];
    subplot(2,2,1)
    histogram(F_noOC,'Normalization','pdf')
    hold on
    [PDFF_noOC,XF_noOC] = ksdensity(F_noOC);
    plot(XF_noOC,PDFF_noOC,'r')
    xlim(xrange)
    ylim(yrange)
    temp = strcat('Factor ', FactorName, ' (No OC)');
    title({'Histogram and Kernel Density Function',temp})

    subplot(2,2,2)
    histogram(F_lev,'Normalization','pdf')
    hold on
    [PDFF_lev,XF_lev] = ksdensity(F_lev);
    plot(XF_lev,PDFF_lev,'r')
    xlim(xrange)
    ylim(yrange)
    temp = strcat('Factor ', FactorName, ' (Lev)');
    title({'Histogram and Kernel Density Function',temp})


    %
    subplot(2,2,3)
    histogram(F_dsg, 'Normalization', 'pdf')
    hold on
    [PDFF_dsg, XF_dsg] = ksdensity(F_dsg);
    plot(XF_dsg, PDFF_dsg, 'r')
    xlim(xrange)
    ylim(yrange)
    temp = strcat('Factor ', FactorName, ' (Dsg)');
    title({'Histogram and Kernel Density Function', temp})
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Inverse Kernel Density Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('**create inverse kernel density functions**\n')
% range
min_prob = 1e-13;
max_prob = 1.0 - min_prob;
num_points = 1e5;
pi = linspace(min_prob, max_prob, num_points);

% No OC
xiNoOC = ksdensity(F_noOC, pi, 'Function', 'icdf');

% lev
xiLev   = ksdensity(F_lev, pi, 'Function', 'icdf');

% dsg
xiDsg   = ksdensity(F_dsg, pi, 'Function', 'icdf');
if do_plots
    figure(2)
    subplot(2,2,1)
    plot(pi, xiNoOC, 'k')
    title({'Inverse CDF distribution No OC', ['Factor ', FactorName, ': ', num2str(N_VP),' Sampled Points']})
    hold on

    subplot(2,2,2)
    plot(pi,xiLev, 'k')
    title({'Inverse CDF distribution Lev', ['Factor ', FactorName,': ',num2str(N_VP),' Sampled Points']})
    hold on

    subplot(2,2,4)
    plot(pi, xiDsg, 'k')
    title({'Inverse CDF distribution Dsg', ['Factor ',FactorName,': ',num2str(N_VP),' Sampled Points']})
    hold on
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Random Paired Samples and Check Desired Diff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LEV
fprintf('create LEV random samples \n')
NUM_TRIALS = 0;
OBJ_flag = 2;

% draw initial samples
xqNoOC          = rand(1,N_VP);
xqLev           = rand(1,N_VP);
samplesNoOC = interp1(pi,xiNoOC,xqNoOC);
samplesLev   = interp1(pi,xiLev,xqLev);

sampleLevdiff  = samplesLev - samplesNoOC;
meanLevdiff      = mean(sampleLevdiff);
stdLevdiff       = std(sampleLevdiff);
OBJ_mean         = abs(meanLevdiff - MEAN_DIFF_LEV);
OBJ_std          = abs(stdLevdiff - SD_DIFF_LEV);



figure;
c = 'red';
ax1 = gca;
plot(ax1, 0,OBJ_std, 'marker', '*', 'color', c)
xlabel(ax1, 'mutation count');
ylabel(ax1,'OBJ std (lev)');
title('LEV standard deviation')
ylim(ax1,[0, OBJ_std]);
hold on

figure;
c2 = 'blue';
ax2 = gca;
plot(ax2, 0,OBJ_mean, 'marker', '*', 'color', c2)
xlabel(ax2, 'mutation count')
ylabel(ax2, 'OBJ mean (lev)')
title('LEV mean')
ylim(ax2, [0, OBJ_mean])
hold on


while( and(NUM_TRIALS < MAX_TRIALS, OBJ_flag > 0))
    NUM_TRIALS = NUM_TRIALS + 1;
    if mod(NUM_TRIALS, 1000) == 0
        fprintf('lev trial number: %i \n', NUM_TRIALS)
        fprintf('obj mean: %d \n',OBJ_mean);
        fprintf('obj std: %d \n',OBJ_std);
    end

    mut_vals = randi(N_VP);
    % take values with biggest difference *from the mean* and redraw
    [~, maxIDs] = maxk(abs((sampleLevdiff-MEAN_DIFF_LEV)), mut_vals);

    new_xqLev         = rand(1,mut_vals);
    new_Lev           = interp1(pi,xiLev,new_xqLev);
    if sum(isnan(new_Lev) > 0)
        fprintf('isnan interpolation, resampling .... \n')
        continue
    end


    new_samplesLev   = samplesLev;

    xqLev_new = xqLev;
    for ii = 1:mut_vals
        vp_val = maxIDs(ii);
        new_samplesLev(vp_val)   = new_Lev(ii);
        xqLev_new(vp_val) = new_xqLev(ii);
    end
    new_sampleLevdiff = new_samplesLev - samplesNoOC;
    %new_meanLevdiff   = mean(new_sampleLevdiff);
    new_stdLevdiff    = std(new_sampleLevdiff);
    %new_OBJ_mean      = abs(new_meanLevdiff - MEAN_DIFF_LEV);
    new_OBJ_std       = abs(new_stdLevdiff - SD_DIFF_LEV);


    % NEW OBJECTIVE: Resample if difference is 2 std away from mean
    resmple_ind = find(abs(new_sampleLevdiff)>(2*SD_DIFF_LEV + abs(MEAN_DIFF_LEV)));

    %fprintf('** Start sample rejection ** \n')
    while length(resmple_ind)>0
        new_xqLev_rej         = rand(1,length(resmple_ind));
        new_Lev_rej           = interp1(pi,xiLev,new_xqLev_rej);
        for jj = 1:length(resmple_ind)
            smpl_ind = resmple_ind(jj);
            new_samplesLev(smpl_ind) = new_Lev_rej(jj);
            xqLev_new(smpl_ind) = new_xqLev_rej(jj);
        end

        new_sampleLevdiff = new_samplesLev - samplesNoOC;

        new_stdLevdiff    = std(new_sampleLevdiff);
        new_OBJ_std       = abs(new_stdLevdiff - SD_DIFF_LEV);



        resmple_ind = find(abs(new_sampleLevdiff)>(2*SD_DIFF_LEV+abs(MEAN_DIFF_LEV)));
    end

    %fprintf('** Finished sample rejection ** \n')




    % keep new samples if OBJ_std improved
    if new_OBJ_std < OBJ_std

        samplesLev   = new_samplesLev;
        xqLev = xqLev_new;


        mut_cnt = mut_cnt + 1;
        %fprintf('mutation %i \n', mut_cnt)

        sampleLevdiff  = samplesLev - samplesNoOC;
        meanLevdiff      = mean(sampleLevdiff);
        stdLevdiff       = std(sampleLevdiff);
        OBJ_mean         = abs(meanLevdiff - MEAN_DIFF_LEV);
        OBJ_std          = abs(stdLevdiff - SD_DIFF_LEV);
        plot(ax1, mut_cnt, OBJ_std, 'marker', '*', 'color', c)
        plot(ax2, mut_cnt, OBJ_mean, 'marker', '*', 'color', c2)
        if mod(mut_cnt, 100) == 0
            pause(2)
        end
    else
        % change per_mutants_lev% of the mutants
        mut_vals = randi(N_VP*per_mutants_lev);
        % take values with biggest difference *from the mean* and redraw
        [~, maxIDs] = maxk(abs((sampleLevdiff-MEAN_DIFF_LEV)), mut_vals);


        new_xqLev         = rand(1,mut_vals);
        new_Lev           = interp1(pi,xiLev,new_xqLev);
        if sum(isnan(new_Lev) > 0)
            fprintf('isnan interpolation, resampling .... \n')
            continue
        end
        new_samplesLev   = samplesLev;

        xqLev_new = xqLev;
        for ii = 1:mut_vals
            vp_vals = maxIDs(ii);
            new_samplesLev(vp_val) = new_Lev(ii);
            xqLev_new(vp_val) = new_xqLev(ii);
        end

        new_sampleLevdiff = new_samplesLev - samplesNoOC;
        new_meanLevdiff   = mean(new_sampleLevdiff);
        %new_stdLevdiff    = std(new_sampleLevdiff);
        new_OBJ_mean      = abs(new_meanLevdiff - MEAN_DIFF_LEV);
        %new_OBJ_std       = abs(new_stdLevdiff - SD_DIFF_LEV);

        if new_OBJ_mean < OBJ_mean
            samplesLev   = new_samplesLev;


            xqLev = xqLev_new;
            mut_cnt = mut_cnt + mut_vals;

            sampleLevdiff  = samplesLev - samplesNoOC;
            meanLevdiff      = mean(sampleLevdiff);
            stdLevdiff       = std(sampleLevdiff);
            OBJ_mean         = abs(meanLevdiff - MEAN_DIFF_LEV);
            OBJ_std          = abs(stdLevdiff - SD_DIFF_LEV);

            plot(ax1, mut_cnt, OBJ_std, 'marker', '*', 'color', c)
            plot(ax2, mut_cnt, OBJ_mean, 'marker', '*', 'color', c2)
            if mod(mut_cnt, 100) == 0
                pause(2)
            end
        end

    end


    sampleLevdiff = samplesLev - samplesNoOC;
    stdLevdiff    = std(sampleLevdiff);
    OBJ_std       = abs(stdLevdiff - SD_DIFF_LEV);


    % NEW OBJECTIVE: Resample if difference is 2 std away from mean
    resmple_ind = find(abs(sampleLevdiff)>(2*SD_DIFF_LEV + abs(MEAN_DIFF_LEV)));

    %fprintf('** Start sample rejection ** \n')
    while length(resmple_ind)>0
        new_xqLev_rej         = rand(1,length(resmple_ind));
        new_Lev_rej           = interp1(pi,xiLev,new_xqLev_rej);
        for jj = 1:length(resmple_ind)
            smpl_ind = resmple_ind(jj);
            samplesLev(smpl_ind) = new_Lev_rej(jj);
            xqLev(smpl_ind) = new_xqLev_rej(jj);
        end

        sampleLevdiff = samplesLev - samplesNoOC;

        stdLevdiff    = std(sampleLevdiff);
        OBJ_std       = abs(stdLevdiff - SD_DIFF_LEV);



        resmple_ind = find(abs(sampleLevdiff)>(2*SD_DIFF_LEV+abs(MEAN_DIFF_LEV)));
    end

    %fprintf('** Finished sample rejection ** \n')
  
    meanLevdiff      = mean(sampleLevdiff);
    stdLevdiff       = std(sampleLevdiff);
    OBJ_mean         = abs(meanLevdiff - MEAN_DIFF_LEV);
    OBJ_std          = abs(stdLevdiff - SD_DIFF_LEV);
    if isnan(OBJ_mean) || isnan(OBJ_std)
        error('value is nan')
    end

    if OBJ_mean <= abs(TOL_lev*MEAN_DIFF_LEV)
        mean_flag = 0;
    else
        mean_flag = 1;
    end
    if OBJ_std <= abs(TOL_lev*SD_DIFF_LEV)
        sd_flag = 0;
    else
        sd_flag = 1;
    end
    OBJ_flag = mean_flag + sd_flag;
end
lev_num_trials = NUM_TRIALS;
fprintf('lev num trials: %i \n', lev_num_trials)

if do_plots
    figure(1)
    subplot(2,2,1)
    histogram(samplesNoOC,'Normalization','pdf')

    subplot(2,2,2)
    histogram(samplesLev,'Normalization','pdf')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DSG
fprintf('\n')
fprintf('create DSG random samples \n')
OBJ_flag  = 2;
NUM_TRIALS = 0;

% draw initial samples
xqDsg           = rand(1,N_VP);
samplesDsg   = interp1(pi,xiDsg,xqDsg);

sampleDsgdiff  = samplesDsg - samplesNoOC;
meanDsgdiff      = mean(sampleDsgdiff);
stdDsgdiff       = std(sampleDsgdiff);
OBJ_mean         = abs(meanDsgdiff - MEAN_DIFF_DSG);
OBJ_std          = abs(stdDsgdiff - SD_DIFF_DSG);

mut_cnt = 0; % mutation count


while( and(NUM_TRIALS < MAX_TRIALS, OBJ_flag > 0))
    NUM_TRIALS = NUM_TRIALS + 1;
    if mod(NUM_TRIALS, 1000) == 0
        fprintf('\t\tdsg trial number: %i \n', NUM_TRIALS)
        fprintf('\t\tobj mean: %d \n',OBJ_mean);
        fprintf('\t\tobj std: %d \n',OBJ_std);
    end
    mut_vals = randi(N_VP);
    % take values with biggest difference *from the mean* and redraw
    [~, maxIDs] = maxk(abs((sampleDsgdiff-MEAN_DIFF_DSG)), mut_vals);

    new_xqDsg             = rand(1,mut_vals);
    new_Dsg           = interp1(pi,xiDsg,new_xqDsg);
    if sum(isnan(new_Dsg) > 0)
        fprintf('isnan interpolation, resampling .... \n')
        continue
    end

    new_samplesDsg   = samplesDsg;
    xqDsg_new = xqDsg;

    for ii = 1:mut_vals
        vp_val = maxIDs(ii);
        new_samplesDsg(vp_val)   = new_Dsg(ii);
        xqDsg_new(vp_val) = new_xqDsg(ii);
    end
    new_sampleDsgdiff = new_samplesDsg - samplesNoOC;
    %new_meanDsgdiff   = mean(new_sampleDsgdiff);
    new_stdDsgdiff    = std(new_sampleDsgdiff);
    %new_OBJ_mean      = abs(new_meanDsgdiff - MEAN_DIFF_DSG);
    new_OBJ_std       = abs(new_stdDsgdiff - SD_DIFF_DSG);

     % NEW OBJECTIVE: Resample if difference is 2 std away from mean
    resmple_ind_dsg = find(abs(new_sampleDsgdiff)>(abs(MEAN_DIFF_DSG) + 2*SD_DIFF_DSG));

 %   fprintf('** Start sample rejection ** \n')
    while length(resmple_ind_dsg)>0
        new_xqDsg_rej         = rand(1,length(resmple_ind_dsg));
        new_Dsg_rej           = interp1(pi,xiDsg,new_xqDsg_rej);
        for jj = 1:length(resmple_ind_dsg)
            smpl_ind = resmple_ind_dsg(jj);
            new_samplesDsg(smpl_ind) = new_Dsg_rej(jj);
            xqDsg_new(smpl_ind) = new_xqDsg_rej(jj);
        end

        new_sampleDsgdiff = new_samplesDsg - samplesNoOC;

        resmple_ind_dsg = find(abs(new_sampleDsgdiff)>(2*SD_DIFF_DSG+abs(MEAN_DIFF_DSG)));


     

        new_stdDsgdiff    = std(new_sampleDsgdiff);
        new_OBJ_std       = abs(new_stdDsgdiff - SD_DIFF_DSG);


    end

   % fprintf('** Finished sample rejection ** \n')


    % keep new samples if OBJ_std improved
    if new_OBJ_std < OBJ_std

        samplesDsg   = new_samplesDsg;

        xqDsg = xqDsg_new;

        mut_cnt = mut_cnt + 1;
        %fprintf('mutation %i \n', mut_cnt)

        sampleDsgdiff  = samplesDsg - samplesNoOC;
        meanDsgdiff      = mean(sampleDsgdiff);
        stdDsgdiff       = std(sampleDsgdiff);
        OBJ_mean         = abs(meanDsgdiff - MEAN_DIFF_DSG);
        OBJ_std          = abs(stdDsgdiff - SD_DIFF_DSG);
        plot(ax1, mut_cnt, OBJ_std, 'marker', '*', 'color', c)
        plot(ax2, mut_cnt, OBJ_mean, 'marker', '*', 'color', c2)
        if mod(mut_cnt, 100) == 0
            pause(2)
        end
    else
        % change 5% of the mutants
        mut_vals = randi(N_VP*per_mutants_dsg);
        % take values with biggest difference *from the mean* and redraw
        [~, maxIDs] = maxk(abs((sampleDsgdiff-MEAN_DIFF_DSG)), mut_vals);
        
        new_xqDsg = rand(1,mut_vals);
        new_Dsg   = interp1(pi,xiDsg,new_xqDsg);
        if sum(isnan(new_Dsg)>0)
            fprintf('isnan interpolation, resampling .... \n')
            continue
        end

        new_samplesDsg = samplesDsg;
        xqDsgnew = xqDsg;

        for ii = 1:mut_vals
            vp_val = maxIDs(ii);
            new_samplesDsg(vp_val) = new_Dsg(ii);
            xqDsg_new(vp_val) = new_xqDsg(ii);
        end


        new_sampleDsgdiff = new_samplesDsg - samplesNoOC;
        new_meanDsgdiff   = mean(new_sampleDsgdiff);
        %new_stdDsgdiff    = std(new_sampleDsgdiff);
        new_OBJ_mean      = abs(new_meanDsgdiff - MEAN_DIFF_DSG);
        %new_OBJ_std       = abs(new_stdDsgdiff - SD_DIFF_DSG);

        if new_OBJ_mean < OBJ_mean
      
            samplesDsg   = new_samplesDsg;
       
            xqDsg = xqDsg_new;


            sampleDsgdiff  = samplesDsg - samplesNoOC;
            meanDsgdiff      = mean(sampleDsgdiff);
            stdDsgdiff       = std(sampleDsgdiff);
            OBJ_mean         = abs(meanDsgdiff - MEAN_DIFF_DSG);
            OBJ_std          = abs(stdDsgdiff - SD_DIFF_DSG);

            plot(ax1, mut_cnt, OBJ_std, 'marker', '*', 'color', c)
            plot(ax2, mut_cnt, OBJ_mean, 'marker', '*', 'color', c2)
            if mod(mut_cnt, 100) == 0
                pause(2)
            end
        end

    end

     sampleDsgdiff = samplesDsg - samplesNoOC;
    stdDsgdiff    = std(sampleDsgdiff);
    OBJ_std       = abs(stdDsgdiff - SD_DIFF_DSG);

     % NEW OBJECTIVE: Resample if difference is 2 std away from mean
    resmple_ind_dsg = find(abs(sampleDsgdiff)>(abs(MEAN_DIFF_DSG) + 2*SD_DIFF_DSG));

 %   fprintf('** Start sample rejection ** \n')
    while length(resmple_ind_dsg)>0
        new_xqDsg_rej         = rand(1,length(resmple_ind_dsg));
        new_Dsg_rej           = interp1(pi,xiDsg,new_xqDsg_rej);
        for jj = 1:length(resmple_ind_dsg)
            smpl_ind = resmple_ind_dsg(jj);
            samplesDsg(smpl_ind) = new_Dsg_rej(jj);
            xqDsg_new(smpl_ind) = new_xqDsg_rej(jj);
        end

        sampleDsgdiff = samplesDsg - samplesNoOC;

        resmple_ind_dsg = find(abs(sampleDsgdiff)>(2*SD_DIFF_DSG+abs(MEAN_DIFF_DSG)));

     

        stdDsgdiff    = std(sampleDsgdiff);
        OBJ_std       = abs(stdDsgdiff - SD_DIFF_DSG);


    end

   % fprintf('** Finished sample rejection ** \n')


            meanDsgdiff      = mean(sampleDsgdiff);
            stdDsgdiff       = std(sampleDsgdiff);
            OBJ_mean         = abs(meanDsgdiff - MEAN_DIFF_DSG);
            OBJ_std          = abs(stdDsgdiff - SD_DIFF_DSG);

    if isnan(OBJ_mean) || isnan(OBJ_std)
        error('value is nan')
    end

    if OBJ_mean <= abs(TOL_lev*MEAN_DIFF_DSG)
        mean_flag = 0;
    else
        mean_flag = 1;
    end
    if OBJ_std <= abs(TOL_lev*SD_DIFF_DSG)
        sd_flag = 0;
    else
        sd_flag = 1;
    end
    OBJ_flag = mean_flag + sd_flag;
end
dsg_num_trials = NUM_TRIALS;
fprintf('dsg num trials: %i \n', dsg_num_trials)


if do_plots

    figure(1)
    subplot(2,2,3)
    histogram(samplesNoOC,'Normalization','pdf')


    subplot(2,2,4)
    histogram(samplesDsg,'Normalization','pdf')
    legend({'data', 'kernel density function', 'virtual patients'})
end




if do_plots

    figure(2)
    subplot(2,2,1)
    plot(xqNoOC,samplesNoOC,'ro')

    subplot(2,2,2)
    plot(xqLev,samplesLev,'ro')
    subplot(2,2,3)
    plot(xqNoOC,samplesNoOC,'bo')

    subplot(2,2,4)
    plot(xqDsg,samplesDsg,'bo')




    figure(7)
    tol_dsg = TOL_dsg;
    tol_lev = TOL_lev;

    NUM_TRIALS_DSG = dsg_num_trials;
    NUM_TRIALS_LEV = lev_num_trials;
    subplot(1,2,2)
    [min_lev,mini_lev] = min(samplesLev-samplesNoOC);
    [max_lev,maxi_lev] = max(samplesLev-samplesNoOC);
    [min_lev,mini_dsg] = min(samplesDsg-samplesNoOC);
    [max_lev,maxi_dsg] = max(samplesDsg-samplesNoOC);


lev_diff = samplesLev-samplesNoOC; 
dsg_diff = samplesDsg-samplesNoOC;

fprintf('Factor %s: Number of DSG samples outside difference range: %d \n',FactorName,sum(abs(dsg_diff) > abs(MEAN_DIFF_DSG)+2*SD_DIFF_DSG)); 
fprintf('Factor %s: Number of LEV samples outside difference range: %d \n',FactorName,sum(abs(lev_diff)>abs(MEAN_DIFF_LEV)+2*SD_DIFF_LEV));



    N = N_VP;
    for i = 1:N
        lineobj = plot([samplesNoOC(i),samplesDsg(i)],'linewidth',1,'color','k');
        hold on
        plot(1, samplesNoOC(i),'.','markersize',15, 'color', lineobj.Color);
        plot(2, samplesDsg(i),'.','markersize',15, 'color', lineobj.Color);
    end

    i = mini_dsg;
    lineobj = plot([samplesNoOC(mini_dsg),samplesDsg(mini_dsg)],'linewidth',1,'color','r');
    plot(1, samplesNoOC(i),'.','markersize',15, 'color', lineobj.Color);
    plot(2, samplesDsg(i),'.','markersize',15, 'color', lineobj.Color);
    i = maxi_dsg;
    lineobj = plot([samplesNoOC(maxi_dsg),samplesDsg(maxi_dsg)],'linewidth',1,'color','g');
    plot(2, samplesDsg(i),'.','markersize',15, 'color', lineobj.Color);
    plot(1, samplesNoOC(i),'.','markersize',15, 'color', lineobj.Color);
    axis([0,3,ylim])
    xticks([1,2])
    xticklabels({'no oc','dsg'})
    set(gcf,'Color','w')
    set(gca,'FontSize',20)
    temp = strcat('Factor ', FactorName, ' %');
    ylabel(temp)
    title({['mean diff: ',num2str(mean(samplesDsg-samplesNoOC))]...
        ,['std diff: ',num2str(std(samplesDsg-samplesNoOC))],['tol: ',num2str(tol_dsg)],['trials: ',num2str(NUM_TRIALS_DSG)]});

    axis([xlim,0,250])
    subplot(1,2,1)
    for i = 1:N
        lineobj = plot([samplesNoOC(i),samplesLev(i)],'linewidth',1,'color','k');
        hold on
        plot(1, samplesNoOC(i),'.','markersize',15, 'color', lineobj.Color);
        plot(2, samplesLev(i),'.','markersize',15, 'color', lineobj.Color);
    end

    i = mini_lev;
    lineobj = plot([samplesNoOC(mini_lev),samplesLev(mini_lev)],'linewidth',1,'color','r');
    plot(1, samplesNoOC(i),'.','markersize',15, 'color', lineobj.Color);
    plot(2, samplesLev(i),'.','markersize',15, 'color', lineobj.Color);
    i = maxi_lev;
    lineobj = plot([samplesNoOC(maxi_lev),samplesLev(maxi_lev)],'linewidth',1,'color','g');
    plot(2, samplesLev(i),'.','markersize',15, 'color', lineobj.Color);
    plot(1, samplesNoOC(i),'.','markersize',15, 'color', lineobj.Color);

    axis([0,3,ylim])
    xticks([1,2])
    xticklabels({'no oc','lev'})
    set(gcf,'Color','w')
    set(gca,'FontSize',20)
    temp = strcat('Factor ', FactorName, ' %')
    ylabel(temp)
    title({['mean diff: ',num2str(mean(samplesLev-samplesNoOC))]...
        ,['std diff: ',num2str(std(samplesLev-samplesNoOC))],['tol: ',num2str(tol_lev)],['trials: ',num2str(NUM_TRIALS_LEV)]});



end
%%%
% save VP
%%%%
if save_VP
    fprintf('saving VP \n')
    samplesMatrix = [samplesNoOC', samplesLev', samplesDsg'];
    T = array2table(samplesMatrix);
    T.Properties.VariableNames = {'noOC', 'lev', 'dsg'};
    FileName = strcat(date, '_Factor',FactorName,'_', 'N_VP-', num2str(N_VP));
    writetable(T, strcat('./samples/', FileName, '.csv'));

    % save information
    savename = strcat('./samples/', FileName, '_info');
    save(savename, 'TOL_lev', 'TOL_dsg', 'lev_num_trials', 'dsg_num_trials', 'MAX_TRIALS',...
        'MEAN_DIFF_LEV', 'MEAN_DIFF_DSG', 'SD_DIFF_LEV', 'SD_DIFF_DSG', 'FactorName')
    fprintf('VP saved to %s \n', FileName)
end
fprintf('simulations complete \n')