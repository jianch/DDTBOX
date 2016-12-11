function [ANALYSIS] = t_tests_classifier_accuracies(ANALYSIS)
%
% This function runs group-level one-sample or paired-samples t tests 
% on classifier accuracy measures derived from each subject. One-sample
% tests are performed against chance level, and paired-samples tests
% compare the decoding accuracy using the observed data with decoding
% accuracy using the permuted data from each subject. The results are
% added to the ANALYSIS structure.
%
% Inputs:
%
%   ANALYSIS            Structure containing organised data and analysis
%                       parameters set in analyse_decoding_erp
%
% Outputs:
%
%   ANALYSIS            Structure containing the data input to the function
%                       plus the results of the statistical analyses.
%
%
% Example:      % [ANALYSIS] = t_tests_classifier_accuracies(ANALYSIS)
%
%
% Copyright (c) 2016 Stefan Bode and contributors
% 
% This file is part of DDTBOX.
%
% DDTBOX is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.



%% Perform group-level tests
for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
        
    for step = 1:size(ANALYSIS.RES.mean_subj_acc,2) % step/time window
            
        % simply test against chance
        if ANALYSIS.permstats == 1
            
            % chance level = 100 / number conditions
            if ANALYSIS.use_robust == 0
                
                [H,P, ~, otherstats] = ttest(ANALYSIS.RES.all_subj_acc(:, na, step), ANALYSIS.chancelevel, ANALYSIS.pstats); % simply against chance
                T = otherstats.tstat;
                clear otherstats;
                
            elseif ANALYSIS.use_robust == 1
                
                % Generate a 'dummy' dataset of chance level values
                chance_level_data_temp = zeros(size(ANALYSIS.RES.all_subj_acc(:, na, step), 1), 1);
                chance_level_data_temp(:) = ANALYSIS.chancelevel;
                % Perform Yuen's t test
                [H,P, ~, T, ~, ~, ~, ~] = yuend_ttest(ANALYSIS.RES.all_subj_acc(:, na, step), chance_level_data_temp, ANALYSIS.trimming, ANALYSIS.pstats);
            
            end % of if ANALYSIS.use_robust
            
        % test against permutation test results    
        elseif ANALYSIS.permstats == 2
            
            % against average permuted distribution
            if ANALYSIS.drawmode == 1

                if ANALYSIS.use_robust == 0
                    
                    [H,P, ~, otherstats] = ttest(ANALYSIS.RES.all_subj_acc(:, na, step), ANALYSIS.RES.all_subj_perm_acc(:,na,step), ANALYSIS.pstats);
                    T = otherstats.tstat;
                    clear otherstats;
                
                elseif ANALYSIS.use_robust == 1
                    
                    [H,P, ~, T, ~, ~, ~, ~] = yuend_ttest(ANALYSIS.RES.all_subj_acc(:, na, step), ANALYSIS.RES.all_subj_perm_acc(:, na, step), ANALYSIS.trimming, ANALYSIS.pstats);
                    
                end % of if ANALYSIS.use_robust
                
            % against one randomly drawn value (from all cross-val repetitions for each participant) for stricter test    
            elseif ANALYSIS.drawmode == 2
                
                for sbj = 1:ANALYSIS.nsbj
                    temp = randperm(size(ANALYSIS.RES.all_subj_perm_acc_reps_draw{sbj,na,step}(:,:),2));
                    drawone = temp(1); clear temp;
                    ANALYSIS.RES.draw_subj_perm_acc(sbj,na,step) = ANALYSIS.RES.all_subj_perm_acc_reps_draw{sbj,na,step}(1,drawone);
                    clear drawone;
                end % sbj
                
                if ANALYSIS.use_robust == 0
                    
                    [H,P, ~, otherstats] = ttest(ANALYSIS.RES.all_subj_acc(:,na,step), ANALYSIS.RES.draw_subj_perm_acc(:,na,step), ANALYSIS.pstats);
                    T = otherstats.tstat;
                    clear otherstats;
                
                elseif ANALYSIS.use_robust == 1
                    
                    [H,P, ~, T, ~, ~, ~, ~] = yuend_ttest(ANALYSIS.RES.all_subj_acc(:,na,step), ANALYSIS.RES.draw_subj_perm_acc(:,na,step), ANALYSIS.trimming, ANALYSIS.pstats);
                    
                end % of if ANALYSIS.use_robust
                
                
            end % if ANALYSIS.drawmode
            
        end % if ANALYSIS.permstats
       
        ANALYSIS.RES.p_ttest(na,step) = P; clear P;
        ANALYSIS.RES.h_ttest_uncorrected(na,step) = H; clear H;
        ANALYSIS.RES.t_ttest(na,step) = T; clear T;
            
    end % of for step = 1:size(ANALYSIS.RES.mean_subj_acc,2) loop
    
end % of for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) loop


%% Correction for multiple comparisons

ANALYSIS.RES.h_ttest = zeros(size(ANALYSIS.RES.mean_subj_acc,1), size(ANALYSIS.RES.mean_subj_acc,2));

switch ANALYSIS.multcompstats

case 0 % No correction for multiple comparisons

    fprintf('\n\nCorrection for multiple comparisons has not been applied\n\n');

    ANALYSIS.RES.h_ttest = ANALYSIS.RES.h_ttest_uncorrected; 

case 1 % Bonferroni Correction

    fprintf('\n\nPerforming corrections for multiple comparisons (Bonferroni)\n\n');

    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
        [ANALYSIS.RES.h_ttest(na, :), ANALYSIS.RES.bonferroni_adjusted_alpha(na)] = multcomp_bonferroni(ANALYSIS.RES.p_ttest(na,:), 'alpha', ANALYSIS.pstats); % Bonferroni correction
        fprintf('The adjusted critical alpha for analysis %i is %1.6f \n', na, ANALYSIS.RES.bonferroni_adjusted_alpha(na));
    end

case 2 % Holm-Bonferroni Correction
    
    fprintf('\n\nPerforming corrections for multiple comparisons (Holm-Bonferroni)\n\n');

    % Here a family of tests is defined as all steps within a given analysis
    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
        [ANALYSIS.RES.h_ttest(na, :), ANALYSIS.RES.holm_adjusted_alpha(na)] = multcomp_holm_bonferroni(ANALYSIS.RES.p_ttest(na,:), 'alpha', ANALYSIS.pstats); % Holm-Bonferroni correction      
        fprintf('\n\nThe adjusted critical alpha for analysis %i is %1.6f   \n\n', na, ANALYSIS.RES.holm_adjusted_alpha(na));
    end % of for na loop    

case 3 % Strong FWER Control Permutation Test (Blaire-Karniski)
    
    fprintf('\n\nPerforming corrections for multiple comparisons (maximum statistic permutation test)\n\n');
    
    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
    
        if ANALYSIS.permstats == 1 % If testing against theoretical chance level
            real_decoding_scores = ANALYSIS.RES.all_subj_acc(:, na, :);
            perm_decoding_scores = zeros(size(real_decoding_scores, 1), 1, size(real_decoding_scores, 3));
            perm_decoding_scores(:, 1, :) = ANALYSIS.chancelevel;
        elseif ANALYSIS.permstats == 2 % If testing against permutation decoding results
            real_decoding_scores = ANALYSIS.RES.all_subj_acc(:, na, :);
            if ANALYSIS.drawmode == 1 % If testing against average permuted distribution
                perm_decoding_scores = ANALYSIS.RES.all_subj_perm_acc(:, na, :);
            elseif ANALYSIS.drawmode == 2 % If testing against one randomly drawn value
                perm_decoding_scores = ANALYSIS.RES.draw_subj_perm_acc(:, na, :);
            end
        end
        
        % Convert to two-dimensional matrix for multcomp correction algorithm
        tmp = squeeze(real_decoding_scores);
        real_decoding_scores = tmp;
        tmp = squeeze(perm_decoding_scores);
        perm_decoding_scores = tmp;

        [ANALYSIS.RES.h_ttest(na, :), ANALYSIS.RES.p_ttest(na,:), ANALYSIS.RES.critical_t(na)] = multcomp_blaire_karniski_permtest(real_decoding_scores, perm_decoding_scores, 'alpha', ANALYSIS.pstats, 'iterations', ANALYSIS.n_iterations);
        fprintf('The adjusted critical t value for analysis %i is %3.3f \n', na, ANALYSIS.RES.critical_t(na));
    end % of for na loop
    
    clear real_decoding_scores
    clear perm_decoding_scores

case 4 % Cluster-Based Permutation Test
 
    fprintf('\n\nPerforming corrections for multiple comparisons (cluster-based permutation test)\n\n');
    
    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
    
        if ANALYSIS.permstats == 1 % If testing against theoretical chance level
            real_decoding_scores = ANALYSIS.RES.all_subj_acc(:, na, :);
            perm_decoding_scores = zeros(size(real_decoding_scores, 1), 1, size(real_decoding_scores, 3));
            perm_decoding_scores(:, 1,:) = ANALYSIS.chancelevel;
        elseif ANALYSIS.permstats == 2 % If testing against permutation decoding results
            real_decoding_scores = ANALYSIS.RES.all_subj_acc(:, na, :);
            if ANALYSIS.drawmode == 1 % If testing against average permuted distribution
                perm_decoding_scores = ANALYSIS.RES.all_subj_perm_acc(:, na, :);
            elseif ANALYSIS.drawmode == 2 % If testing against one randomly drawn value
                perm_decoding_scores = ANALYSIS.RES.draw_subj_perm_acc(:, na, :);
            end    
        end
        % Convert to two-dimensional matrix for multcomp correction algorithm
        tmp = squeeze(real_decoding_scores);
        real_decoding_scores = tmp;
        tmp = squeeze(perm_decoding_scores);
        perm_decoding_scores = tmp;
        
        [ANALYSIS.RES.h_ttest(na, :)] = multcomp_cluster_permtest(real_decoding_scores, perm_decoding_scores,  'alpha', ANALYSIS.pstats, 'iterations', ANALYSIS.n_iterations, 'clusteringalpha', ANALYSIS.cluster_test_alpha);
    end % of for na loop
    clear real_decoding_scores
    clear perm_decoding_scores    

case 5 % KTMS Generalised FWER Control Using Permutation Testing
    
    fprintf('\n\nPerforming corrections for multiple comparisons (KTMS generalised FWER control)\n\n');

    % Adapted from permutation test script
    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis

        if ANALYSIS.permstats == 1 % If testing against theoretical chance level
            real_decoding_scores = ANALYSIS.RES.all_subj_acc(:, na, :);
            perm_decoding_scores = zeros(size(real_decoding_scores, 1), 1, size(real_decoding_scores, 3));
            perm_decoding_scores(:, 1,:) = ANALYSIS.chancelevel;
        elseif ANALYSIS.permstats == 2 % If testing against permutation decoding results
            real_decoding_scores = ANALYSIS.RES.all_subj_acc(:, na, :);
            if ANALYSIS.drawmode == 1 % If testing against average permuted distribution
                perm_decoding_scores = ANALYSIS.RES.all_subj_perm_acc(:, na, :);
            elseif ANALYSIS.drawmode == 2 % If testing against one randomly drawn value
                perm_decoding_scores = ANALYSIS.RES.draw_subj_perm_acc(:, na, :);
            end
        end
        % Convert to two-dimensional matrix for multcomp correction algorithm
        tmp = squeeze(real_decoding_scores);
        real_decoding_scores = tmp;
        tmp = squeeze(perm_decoding_scores);
        perm_decoding_scores = tmp;

        [ANALYSIS.RES.h_ttest(na, :)] = multcomp_ktms(real_decoding_scores, perm_decoding_scores, 'alpha', ANALYSIS.pstats, 'iterations', ANALYSIS.n_iterations, 'ktms_u', ANALYSIS.ktms_u);
    end % of for na loop
    clear real_decoding_scores
    clear perm_decoding_scores

case 6 % Benjamini-Hochberg FDR Control

    fprintf('\n\nPerforming corrections for multiple comparisons (Benjamini-Hochberg)\n\n');
    
    % Here a family of tests is defined as all steps within a given analysis
    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
        [ANALYSIS.RES.h_ttest(na, :), ANALYSIS.RES.bh_crit_alpha(na)] = multcomp_fdr_bh(ANALYSIS.RES.p_ttest(na,:), 'alpha', ANALYSIS.pstats);
        fprintf('\n\nThe adjusted critical alpha for analysis %i is %1.6f \n\n', na, ANALYSIS.RES.bh_crit_alpha(na));
    end % of for na loop

case 7 % Benjamini-Krieger-Yekutieli FDR Control
    
    fprintf('\n\nPerforming corrections for multiple comparisons (Benjamini-Krieger-Yekutieli)\n\n');

    % Here a family of tests is defined as all steps within a given analysis
    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
        [ANALYSIS.RES.h_ttest(na, :), ANALYSIS.RES.bky_crit_alpha(na)] = multcomp_fdr_bky(ANALYSIS.RES.p_ttest(na,:), 'alpha', ANALYSIS.pstats);
        fprintf('The adjusted critical alpha for analysis %i is %1.6f \n', na, ANALYSIS.RES.bky_crit_alpha(na));
    end % of for na loop

case 8 % Benjamini-Yekutieli FDR Control
    
    fprintf('\n\nPerforming corrections for multiple comparisons (Benjamini-Yekutieli)\n\n');

    % Here a family of tests is defined as all steps within a given analysis
    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
        [ANALYSIS.RES.h_ttest(na, :), ANALYSIS.RES.by_crit_alpha(na)] = multcomp_fdr_by(ANALYSIS.RES.p_ttest(na,:), 'alpha', ANALYSIS.pstats);
        fprintf('The adjusted critical alpha for analysis %i is %1.6f \n', na, ANALYSIS.RES.by_crit_alpha(na));
    end % of for na loop

% If some other option is chosen then do not correct for multiple comparisons, but notify user
otherwise
    fprintf('\n\nUnavailable multiple comparisons option chosen. Will use uncorrected p-values \n\n');
    ANALYSIS.RES.h_ttest = ANALYSIS.RES.h_ttest_uncorrected; 
end % of ANALYSIS.multcompstats switch

% Marking h values (statistical significance) for plotting
ANALYSIS.RES.h = ANALYSIS.RES.h_ttest;
