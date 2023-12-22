%Clear
clear all
close all

% Goal:
% Create a sample of individual with Factor II levels
% that come from the observed data for individual factors 
% Middeldorp et al. 2000 (Figure 2) and the difference
% in mean reported in Table 1

% Goal: There is a SHIFT for each patient. The patient under treatment
% is guaranteed to remain above or below the mean effect. 

% Dec 19, 2023
% Problem: VP factor level changes were fully independent between OCs.
%          This meant that a patients thrombin response was uncorrelated
%          between OCs. (Example: They could be Higher in Factor II in OC1
%          and far lower in Factor II in OC2).
% Idea:    (1) Patients in the lower (or upper) half before OC chould
%          be restricted to stay in that lower (or upper) half
%          (2) This range of values (upper or lower) chould be same for 
%          factor levels between OCs. 
% Coding Idea: We group patients into up or down response for a factor
%          and then relate that response across OCs.
%          Pick under OC1 the response (independently for each patient)
%          Assume the direction of change will be same under OC2.
% Even Better Idea:
%          Hyperparamters: Patients in the top 1% of a response on OC1
%          would also be in the top 1% of a response in OC2.
%          Hyperparameter: Bias in the direction of change. Right now
%          unbiased. We could say if you're at the top 10% before OC you
%          should be within the top 20% after OC. 
%          Potentially look at the Kendal-Tau distance
%          (https://en.wikipedia.org/wiki/Kendall_tau_distance).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data and Set Ranges %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load Data
csvtable1 =  readtable('data/FactorII_noOC1.csv');
csvtable2 =  readtable('data/FactorII_lev.csv');
FactorII_noOC1 = csvtable1.Var2; 
FactorII_lev   = csvtable2.Var2;

mean(FactorII_noOC1);
std(FactorII_noOC1);

mean(FactorII_lev);
std(FactorII_lev);

%Desired Mean Difference After Treatment (Lev - NoOC):
MEAN = 12; %Table 1 taken from Middeldorp et al. 2000
STD  = 8;  %The Standard Deviation.

%How many patients:
N = 1000; 

%How many trials do we want:
MAX_TRIALS = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Plots of Original Distributions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Plot Kernel Density Function
% subplot(1,2,1)
% histogram(FactorII_noOC1,'Normalization','pdf')
% hold on
% [PDFFactorII_noOC1,XFactorII_noOC1] = ksdensity(FactorII_noOC1);
% plot(XFactorII_noOC1,PDFFactorII_noOC1,'r')
% title({'Histogram and Kernel Density Function','Factor II (No OC)'})
% 
% subplot(1,2,2)
% histogram(FactorII_lev,'Normalization','pdf')
% hold on
% [PDFFactorII_lev,XFactorII_lev] = ksdensity(FactorII_lev);
% plot(XFactorII_lev,PDFFactorII_lev,'r')
% title({'Histogram and Kernel Density Function','Factor II (Lev)'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Inverse Kernel Density Functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Range 
pi = linspace(.0001,.9999,10000);

%No OC Factor: 
xiNoOC1 = ksdensity(FactorII_noOC1,pi,'Function','icdf');

%Lev Factors
xiLev = ksdensity(FactorII_lev,pi,'Function','icdf');

% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Random Paired Samples and Check Desired Diff %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OBJ        = 0; %Do we match the objective. 
NUM_TRIALS = 0;

%Pick your NO OC Values for Each Patient
xqNoOC1     = rand(1,N);
samplesNoOC = interp1(pi,xiNoOC1,xqNoOC1);

%Min Shift:
meanNoOC = mean(FactorII_noOC1);

minRange = samplesNoOC + MEAN - (samplesNoOC<meanNoOC) .*(3*STD);%*rand(1,N));
maxRange = samplesNoOC + MEAN + (samplesNoOC>meanNoOC) .*(3*STD);%*rang(1,N));

%If a patient is BELOW the mean, 
%(*) min value will be noOC + MEAN - 2*STD
%(*) max value will be noOC + MEAN 

%If a patient is ABOVE the mean,
%(*) min value will be noOC + MEAN ;
%(*) max value will be noOC + MEAN + 2*STD

%Determine the min/max value under OC for each patient
%Assume we stay 2 standard deviations away from the mean.
%minRange = samplesNoOC - (MEAN + STD*rand(1,N)); %Should check we stay ABOVE minLev;
%maxRange = samplesNoOC + (MEAN + STD*rand(1,N)); %Must check we stay BELOW the maxLev;

%Check that the min/max for each individual is within the
%domain for the kernel density functions for the lev distribution
for i=1:N
    if(minRange(i)<min(xiLev))
        minRange(i) = min(xiLev)*1.01;
    end
    if(maxRange(i)>max(xiLev))
        maxRange(i) = max(xiLev)*.99;
    end
end

%Translate these min and max values to the percentiles of the lev KDF
minPct   =  interp1(xiLev,pi,minRange);
maxPct   =  interp1(xiLev,pi,maxRange);

%Generate the appropriate range for the percentile. 
xqLev       = minPct + (maxPct - minPct).*rand(1,N); 
samplesLev  = interp1(pi,xiLev,xqLev);


% while ( and(OBJ > 1, NUM_TRIALS < MAX_TRIALS ) ) 
% 
%     NUM_TRIALS = NUM_TRIALS + 1
% 
%     %(a) Create two LHS Samples in [0,1] (keeps people in same order)
%     %xqNoOC1 = lhsdesign(1,N); 
%     %xqLev   = lhsdesign(1,N);
% 
%     %(b) Create a MonteCarlo Sample in [minPct,maxPct] (Each Uniform in [0,1])
%     %xqLev   = rand(1,N);
%     xqLev    = minPct + (maxPct - minPct).*rand(1,N); 
%     a + (b-a).*rand(100,1)
%     %Convert each of the [0,1] into [min,max]
% 
%     samplesLev  = interp1(pi,xiLev,xqLev);
% 
%     meanSampleDiff = mean(samplesLev - samplesNoOC);
%     stdSampleDiff  = std(samplesLev - samplesNoOC);
%     
%     OBJ = abs(MEAN_DIFF - meanSampleDiff);
%     
% end
% 
% % Alternative Ideas
% % (1)Genetic Algorithm;
% %   * Generate a random *New Sample* & Pick an Existing Sample 
% %   * If the *New Sample* is closer MEAN_DIFF than the Existing Sample, replace, 
% % (2) Permutation
% %   * We know the data was from the same individuals, we do not
% %     know which data points should be paired. Try it.
% 
% 
% % OBJ        = 2*RANGE;
% % NUM_TRIALS = 0;
% % 
% % while ( and(OBJ > RANGE, NUM_TRIALS < MAX_TRIALS ) ) 
% %     
% %     NUM_TRIALS = NUM_TRIALS + 1
% %     
% %     %Randomly Permute Patients:
% %     perm = randperm(28)
% % 
% %     %Pick Paired Samples
% %     samplesNoOC = FactorII_noOC1;
% %     samplesLev  = FactorII_lev(perm);
% % 
% %     meanSampleDiff = mean(samplesLev - samplesNoOC);
% %     
% %     OBJ = abs(MEAN_DIFF - meanSampleDiff);
% %     
% % end
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot Final Sets of Samples %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure(2)
% subplot(1,2,1)
% plot(xqNoOC1,samplesNoOC,'ro')
% 
% subplot(1,2,2)
% plot(xqLev,samplesLev,'ro')
% 
% %figure(3)
% figure(1)
% subplot(1,2,1)
% histogram(samplesNoOC,'Normalization','pdf')
% title('Histogram of Sampled Factor II No OC')
% 
% subplot(1,2,2)
% histogram(samplesLev,'Normalization','pdf')
% title('Histogram of Sampled Factor II Lev')
