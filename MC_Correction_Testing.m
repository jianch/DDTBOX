% MC_Correction_Testing.m
%
% Implements and tests multiple comparisons corrections before implementing
% these in DDTBOX. We can also run simulations to ensure that each method
% controls the FWER or FDR at the nominal level.
% Written by Daniel Feuerriegel on 4/16 as part of the development of DDTBOX.
%


%% Housekeeping
clear all;
close all;


%% Setting variables (number of simulation subjects etc.)
nTests = 40; % Number of tests to run
nSubjects = 40; 
meanDiff = 0; % Mean effect magnitude
SD = 1; % SD of Gaussian distribution

alphaLevel = 0.05; % Critical p-value

nIterations = 100; % Number of bootstrap samples to take in strong FWER control permutation testing





%% Make sets of data (each as "step" in DDTBox).
sample = zeros(nTests, nSubjects);
for test = 1:nTests
    sample(test, 1:nSubjects) = SD * randn(1, nSubjects) + meanDiff;
    
    [h(test), p(test)] = ttest(sample(test,1:nSubjects));
    
end

%% Bonferroni correction
bonferroni_corrected_alpha = alphaLevel / nTests;
for test = 1:nTests
    if p(test) < bonferroni_corrected_alpha   
        bonferroni_h(test) = 1;
    else
        bonferroni_h(test) = 0;
    end
end


%% Holm-Bonferroni Correction
% Sort p-values from smallest to largest
sorted_p = sort(p);
foundCritAlpha = 0; % Reset to signify we have not found critical alpha level
% Find critical alpha level (sig. for all smaller p-values than this)
for holmStep = 1:nTests;
   if sorted_p(holmStep) > alphaLevel / (nTests + 1 - holmStep) & foundCritAlpha == 0
       holm_corrected_alpha = sorted_p(holmStep);
       foundCritAlpha = 1;
   end
end

if ~exist('holm_corrected_alpha') % If all null hypotheses are rejected
    holm_corrected_alpha = alphaLevel;
end
% Declare tests significant if they are smaller than the adjusted critical alpha
for test = 1:nTests
    if p(test) < holm_corrected_alpha   
        holm_h(test) = 1;
    else
        holm_h(test) = 0;
    end
end


%% Benjamini-Hochberg Procedure

% Sort p-values from smallest to largest
sorted_p = sort(p);
% Find critical k value (see tutorial notes)
for BenHoch_Step = 1:nTests;
    if sorted_p(BenHoch_Step) <= (BenHoch_Step / nTests) * alphaLevel
        BenHoch_critical_alpha = sorted_p(BenHoch_Step);
    end
end
% If no tests are significant set critical alpha to zero
if ~exist('BenHoch_critical_alpha')
    BenHoch_critical_alpha = 0;
end

% Declare tests significant if they are smaller than or equal to the adjusted critical alpha
for test = 1:nTests
    if p(test) <= BenHoch_critical_alpha   
        BenHoch_h(test) = 1;
    else
        BenHoch_h(test) = 0;
    end
end


%% Benjamini-Yekutieli Procedure

% Sort p-values from smallest to largest
sorted_p = sort(p);

% To help calculate the Benjamini-Yekutieli critical alpha
jValues = zeros(1,nTests);
for j = 1:nTests
    jValues(j) = 1 / j;
end
% Find critical k value (see tutorial notes)
for BenYek_Step = 1:nTests;
    if sorted_p(BenYek_Step) <= (BenYek_Step / nTests * sum(jValues)) * alphaLevel
        BenYek_critical_alpha = sorted_p(BenYek_Step);
    end
end
% If no tests are significant set critical alpha to zero
if ~exist('BenYek_critical_alpha')
    BenYek_critical_alpha = 0;
end

% Declare tests significant if they are smaller than or equal to the adjusted critical alpha
for test = 1:nTests
    if p(test) <= BenYek_critical_alpha   
        BenYek_h(test) = 1;
    else
        BenYek_h(test) = 0;
    end
end



%% Benjamini-Krieger-Yekutieli Procedure
% Stage 1: Estimate the number of false null hypotheses using the modified
% alpha level.

% Sort p-values from smallest to largest
sorted_p = sort(p);
% Find critical k value (see tutorial notes)
for BKY_Step = 1:nTests;
    if sorted_p(BKY_Step) <= (BKY_Step / nTests) * (alphaLevel / ( 1 + alphaLevel));
        BKY_Stage1_critical_alpha = sorted_p(BKY_Step);
    end
end
% If no tests are significant set critical alpha to zero
if ~exist('BKY_Stage1_critical_alpha')
    BKY_Stage1_critical_alpha = 0;
end

% Declare tests significant if they are smaller than or equal to the adjusted critical alpha
for test = 1:nTests
    if p(test) <= BKY_Stage1_critical_alpha   
        BKY_Stage1_h(test) = 1;
    else
        BKY_Stage1_h(test) = 0;
    end
end

% Count the number of rejected null hypotheses (for step 2)
BKY_Stage1_nRejections = sum(BKY_Stage1_h);

if BKY_Stage1_nRejections == 0; % if no null hypotheses were rejected
    BKY_Stage2_h(1:nTests) = 0; % Don't reject any hypotheses
elseif BKY_Stage1_nRejections == nTests; % if all null hypotheseses were rejected
    BKY_Stage2_h(1:nTests) = 1; % Reject all hypotheses
else % If some (but not all) null hypotheses were rejected  
    
    for BKY_Step = 1:nTests;
        if sorted_p(BKY_Step) <= (BKY_Step / nTests) * ( (nTests / (nTests - BKY_Stage1_nRejections) ) * (alphaLevel / ( 1 + alphaLevel)) );
            BKY_Stage2_critical_alpha = sorted_p(BKY_Step);
        end
    end

% If no tests are significant set critical alpha to zero
if ~exist('BKY_Stage2_critical_alpha')
    BKY_Stage2_critical_alpha = 0;
end

% Declare tests significant if they are smaller than or equal to the adjusted critical alpha
for test = 1:nTests
    if p(test) <= BKY_Stage2_critical_alpha   
        BKY_Stage2_h(test) = 1;
    else
        BKY_Stage2_h(test) = 0;
    end
end

end


%% Strong FWER Control Permutation Testing

% NOTE: Permutation testing as outlined in Groppe et al. (2011) is not
% suitable for the DDTBox corrections, as the permutation testing in Groppe
% et al. works with label switching between two conditions (permutation
% testing at the group level), whereas in DDTBOX permutation test decoding is done 
% at the subject level. For this reason I will try using bootstrapping of
% the subject-level permutation test results instead.
%
% TODO: Check Efron & Tibshirani book to see how they deal with one-sample
% tests.

% TODO: Generate 'real' and 'null' samples at the start of the script for
% use here.


% Generate t(max) distribution from the null data (use permutation results
% when implementing in DDTBOX)
t_max = zeros(1, nIterations);

for iteration = 1:nIterations
    
    % Draw a random sample for each test
    for test = 1:nTests
        temp(test, 1:nSubjects) = randsample(sample(test, 1:nSubjects), nSubjects, true); % Draw a bootstrap sample (i.e. with replacement)
        [~, ~, ~, temp_stats] = ttest(temp(test, 1:nSubjects));
        t_stat(test, iteration) = temp_stats.tstat;
    end    
    
    % Get the maximum t-value within the family of tests and store in a
    % vector. This is to create a null hypothesis distribution.
    t_max(iteration) = max(abs(t_stat(:, iteration)));
end










%% Calculate the FWER/FDR of the tests
FWER.uncorrected = sum(h) / nTests;
FWER.bonferroni = sum(bonferroni_h) / nTests;
FWER.holm = sum(holm_h) / nTests;
FWER.BenHoch = sum(BenHoch_h) / nTests;
FWER.BenYek = sum(BenYek_h) / nTests;
FWER.BKY = sum(BKY_Stage2_h) / nTests;