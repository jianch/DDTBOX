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
nSubjects = 20; 
meanDiff = 0; % Mean effect magnitude
SD = 1; % SD of Gaussian distribution

alphaLevel = 0.05; % Critical p-value

nIterations = 100; % Number of bootstrap samples to take in strong FWER control permutation testing

KTMS_u = 1; % u parameter of the KTMS GFWER control procedure

% Extra distributions to make for permutation-based procedures
nRealEffects = 8; % Number of real effects in the data.
meanDiffRealEffects = 1; % Mean effect magnitude



%% Make sets of data (each as "step" in DDTBox).
sample = zeros(nTests, nSubjects);
for test = 1:nTests
    sample(test, 1:nSubjects) = SD * randn(1, nSubjects) + meanDiff;
    
    [h(test), p(test)] = ttest(sample(test,1:nSubjects));
    
end

% Make data showing real effects to use in permutation tests
sampleForRealEffects = zeros(nTests, nSubjects); % Preallocate matrix
realEffectTestIndices = randi(nTests, 1, nRealEffects); % Randomly allocate location of real effects
% Mark locations of the true effects
realEffectLocations = zeros(1, nTests);
realEffectLocations(realEffectTestIndices) = 1;

% Generate null samples in the same way as above
for test = 1:nTests
    sampleForRealEffects(test, 1:nSubjects) = SD * randn(1, nSubjects) + meanDiff;    
end
% Add a fixed value to a subset of tests (adding some real effects into the data)
sampleForRealEffects(realEffectTestIndices, :) =  sampleForRealEffects(realEffectTestIndices, :) + meanDiffRealEffects;

for test = 1:nTests    
    [realEffect_h(test), realEffect_p(test), ~, tempStats] = ttest(sampleForRealEffects(test, 1:nSubjects), 0, 'Alpha', alphaLevel);
    realEffect_t(test) = tempStats.tstat;
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
% NOTE: The FWER control of this method has not been validated as far as I
% know (although it may be somewhere in the literature). Therefore we will
% want to run simulations to ensure that it does control the FWER.
%
%
% TODO: Check Efron & Tibshirani book to see how they deal with one-sample
% tests.


% Generate t(max) distribution from the null data (use permutation results
% when implementing in DDTBOX)
t_max = zeros(1, nIterations);

for iteration = 1:nIterations
    
    % Draw a random sample for each test
    for test = 1:nTests
        temp(test, 1:nSubjects) = randsample(sample(test, 1:nSubjects), nSubjects, true); % Draw a bootstrap sample (i.e. with replacement)
        [~, ~, ~, temp_stats] = ttest(temp(test, 1:nSubjects), 0, 'Alpha', alphaLevel);
        t_stat(test, iteration) = abs(temp_stats.tstat);
    end    
    
    % Get the maximum t-value within the family of tests and store in a
    % vector. This is to create a null hypothesis distribution.
    t_max(iteration) = max(t_stat(:, iteration));
end

% Calculating the 95th percentile of t_max values (used as decision
% critieria for statistical significance)
PermTest_Null_Cutoff = prctile(t_max, ((1 - alphaLevel) * 100));

% Checking whether each test statistic is above the specified threshold:
for test = 1:nTests
    if abs(realEffect_t(test)) > PermTest_Null_Cutoff;
        permTest_h(test) = 1;
    else
        permTest_h(test) = 0;
    end
end



%% Cluster-Based Permutation Testing (Weak FWER Control)

% Generate a maximum cluster mass distribution from the null data (use permutation results
% when implementing in DDTBOX)
max_cluster_mass = zeros(1, nIterations);
t_stat_for_clustering = zeros(nTests, nIterations);
clusterPermTest_h = zeros(nTests, nIterations);

for iteration = 1:nIterations
    
    % Draw a random bootstrap sample for each test
    for test = 1:nTests
        temp(test, 1:nSubjects) = randsample(sample(test, 1:nSubjects), nSubjects, true); % Draw a bootstrap sample (i.e. with replacement)
        [clusterPermTest_h(test, iteration), ~, ~, temp_stats] = ttest(temp(test, 1:nSubjects), 0, 'Alpha', alphaLevel);
        t_stat_for_clustering(test, iteration) = temp_stats.tstat;
    end    
    
    % Identify clusters and generate a maximum cluster statistic
    clusterMassVector = [0]; % Resets vector of cluster masses
    clusterCounter = 0;
    for test = 1:nTests     
        if clusterPermTest_h(test, iteration) == 1
            if test == 1 % If the first test in the set
                clusterCounter = clusterCounter + 1;
                clusterMassVector(clusterCounter) = abs(t_stat_for_clustering(test, iteration));
            else
                if clusterPermTest_h(test - 1, iteration) == 1
                    clusterMassVector(clusterCounter) = clusterMassVector(clusterCounter) + abs(t_stat_for_clustering(test, iteration));
                elseif clusterPermTest_h(test - 1, iteration) == 0
                    clusterCounter = clusterCounter + 1;
                    clusterMassVector(clusterCounter) = t_stat_for_clustering(test, iteration);
                end 
            end % of if test == 1
        end % of if clusterPermTest
    end % of for tests = 1:nTests
    
    % Find the maximum cluster mass
    max_cluster_mass(iteration) = max(clusterMassVector);
end

% Calculating the 95th percentile of maximum cluster mass values (used as decision
% critieria for statistical significance)
clusterMass_Null_Cutoff = prctile(max_cluster_mass, ((1 - alphaLevel) * 100));

% Calculate cluster masses in the actual (non-permutation) tests
clusterMassVector = [0]; % Resets vector of cluster masses
clusterCounter = 0;
clusterLocations = zeros(1, nTests);
clusterCorrected_Sig_Tests = zeros(1, nTests);

for test = 1:nTests    
    [realEffect_h(test), realEffect_p(test), ~, tempStats] = ttest(sampleForRealEffects(test, 1:nSubjects), 0, 'Alpha', alphaLevel);
    realEffect_t(test) = tempStats.tstat;
    
    if realEffect_h(test) == 1
        if test == 1 % If the first test in the set
            clusterCounter = clusterCounter + 1;
            clusterMassVector(clusterCounter) = abs(realEffect_t(test));
            clusterLocations(test) = clusterCounter;
            % Tagging as positive or negative sign effect
            if realEffect_t < 0
                t_sign(test) = -1;
            else
                t_sign(test) = 1;
            end
            
        elseif test > 1
            % Tagging as positive or negative sign effect
            if realEffect_t < 0
                t_sign(test) = -1;
            else
                t_sign(test) = 1;
            end
            % Add to the same cluster only if the previous test was sig.
            % and of the same sign (direction).
            if realEffect_h(test - 1) == 1 && t_sign(test - 1) == t_sign(test)
                    clusterMassVector(clusterCounter) = clusterMassVector(clusterCounter) + abs(realEffect_t(test));
                    clusterLocations(test) = clusterCounter;
            else
                clusterCounter = clusterCounter + 1;
                clusterMassVector(clusterCounter) = abs(realEffect_t(test));
                clusterLocations(test) = clusterCounter;
            end 
        end % of if test == 1
    end % of if clusterPermTest  
end % of for test = 1:nTests

for clusterNo = 1:length(clusterMassVector);
    if clusterMassVector(clusterNo) > clusterMass_Null_Cutoff
        clusterCorrected_Sig_Tests(clusterLocations == clusterNo) = 1;
    end
end


%% KTMS Generalised FWER Control Using Permutation Testing

KTMS_sig_effect_locations = zeros(1, nTests);
% Sort p-values from smallest to largest
sorted_p = sort(realEffect_p);
% Automatically reject the u smallest hypotheses (u is set by user as KTMS_u variable).
KTMS_autoRejectAlpha = sorted_p(KTMS_u);

KTMS_sig_effect_locations(realEffect_p <= KTMS_autoRejectAlpha) = 1; % Mark tests with u smallest p-values as statistically significant.

% Run strong FWER control permutation test but use u + 1th most extreme
% test statistic.
KTMS_t_max = zeros(1, nIterations);
t_stat = zeros(nTests, nIterations);

for iteration = 1:nIterations
    
    % Draw a random sample for each test
    for test = 1:nTests
        temp(test, 1:nSubjects) = randsample(sample(test, 1:nSubjects), nSubjects, true); % Draw a bootstrap sample (i.e. with replacement)
        [~, ~, ~, temp_stats] = ttest(temp(test, 1:nSubjects), 0, 'Alpha', alphaLevel);
        t_stat(test, iteration) = abs(temp_stats.tstat);
    end    
    
    % Get the maximum t-value within the family of tests and store in a
    % vector. This is to create a null hypothesis distribution.
    t_sorted = sort(t_stat(:, iteration), 'descend');
    
    KTMS_t_max(iteration) = t_sorted(KTMS_u + 1);
end

% Calculating the 95th percentile of t_max values (used as decision
% critieria for statistical significance)
KTMS_Null_Cutoff = prctile(KTMS_t_max, ((1 - alphaLevel) * 100));


% Checking whether each test statistic is above the specified threshold:
for test = 1:nTests
    if abs(realEffect_t(test)) > KTMS_Null_Cutoff;
        KTMS_sig_effect_locations(test) = 1;
    end
end









%% Comparing Permutation Test Outcomes
% Show tests with real effects and the results of t(max) and
% cluster-corrected permutatation tests:
figure;
A(1,1:nTests) = realEffectLocations;
A(2,1:nTests) = permTest_h;
A(3,1:nTests) = clusterCorrected_Sig_Tests;
A(4,1:nTests) = KTMS_sig_effect_locations;
imagesc(A);



%% Calculate the FWER/FDR of the tests
FWER.uncorrected = sum(h) / nTests;
FWER.bonferroni = sum(bonferroni_h) / nTests;
FWER.holm = sum(holm_h) / nTests;
FWER.BenHoch = sum(BenHoch_h) / nTests;
FWER.BenYek = sum(BenYek_h) / nTests;
FWER.BKY = sum(BKY_Stage2_h) / nTests;

% Calculate how many tests were found to be statistically-significant with
% mixed true + false effects
NSigTests.permTest = sum(permTest_h) / nTests; 
NSigTests.clusterPermTest = sum(clusterCorrected_Sig_Tests) / nTests; 
NSigTests.KTMS_sig_effect_locations = sum(KTMS_sig_effect_locations) / nTests; 