
function [~] = t_tests_classifier_accuracies(input1, input2, etc)
% Put header here
%
%

for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
        
    for step = 1:size(ANALYSIS.RES.mean_subj_acc,2) % step
            
        % simply test against chance
        if ANALYSIS.permstats == 1
            
            % chance level = 100 / number conditions
            [H,P, ~, otherstats] = ttest(ANALYSIS.RES.all_subj_acc(:,na,step),ANALYSIS.chancelevel,ANALYSIS.pstats); % simply against chance
            T = otherstats.tstat;
            clear otherstats;
            
        % test against permutation test results    
        elseif ANALYSIS.permstats == 2
            
            % against average permuted distribution
            if ANALYSIS.drawmode == 1
                
                [H,P, ~, otherstats] = ttest(ANALYSIS.RES.all_subj_acc(:,na,step),ANALYSIS.RES.all_subj_perm_acc(:,na,step),ANALYSIS.pstats);
                T = otherstats.tstat;
                clear otherstats;
                
            % against one randomly drawn value (from all cross-val repetitions for each participant) for stricter test    
            elseif ANALYSIS.drawmode == 2
                
                for sbj = 1:ANALYSIS.nsbj
                    temp = randperm(size(ANALYSIS.RES.all_subj_perm_acc_reps_draw{sbj,na,step}(:,:),2));
                    drawone = temp(1); clear temp;
                    ANALYSIS.RES.draw_subj_perm_acc(sbj,na,step) = ANALYSIS.RES.all_subj_perm_acc_reps_draw{sbj,na,step}(1,drawone);
                    clear drawone;
                end % sbj
                
                [H,P, ~, otherstats] = ttest(ANALYSIS.RES.all_subj_acc(:,na,step),ANALYSIS.RES.draw_subj_perm_acc(:,na,step),ANALYSIS.pstats);
                T = otherstats.tstat;
                clear otherstats;
            end % if ANALYSIS.drawmode
            
        end % if ANALYSIS.permstats
       
        ANALYSIS.RES.p_ttest(na,step) = P; clear P;
        ANALYSIS.RES.h_ttest_uncorrected(na,step) = H; clear H;
        ANALYSIS.RES.t_ttest(na,step) = T; clear T;
            
    end % of for step = 1:size(ANALYSIS.RES.mean_subj_acc,2) loop
    
end % of for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) loop


%% CORRECTION FOR MULTIPLE COMPARISONS
%__________________________________________________________________________

ANALYSIS.RES.h_ttest = zeros(size(ANALYSIS.RES.mean_subj_acc,1), size(ANALYSIS.RES.mean_subj_acc,2));

switch ANALYSIS.multcompstats

case 0 % No correction for multiple comparisons

    fprintf('Correction for multiple comparisons has not been applied\n');

    ANALYSIS.RES.h_ttest = ANALYSIS.RES.h_ttest_uncorrected; 

%__________________________________________________________________________

case 1 % Bonferroni Correction

    fprintf('Performing corrections for multiple comparisons (Bonferroni)\n');

    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
        [ANALYSIS.RES.h_ttest(na, :), ANALYSIS.RES.bonferroni_adjusted_alpha(na)] = multcomp_bonferroni(ANALYSIS.RES.p_ttest(na,:), 'alpha', ANALYSIS.pstats); % Bonferroni correction
        fprintf('The adjusted critical alpha for analysis %i is %1.6f \n', na, ANALYSIS.RES.bonferroni_adjusted_alpha(na));
    end

%__________________________________________________________________________    

case 2 % Holm-Bonferroni Correction
    
    fprintf('Performing corrections for multiple comparisons (Holm-Bonferroni)\n');

    % Here a family of tests is defined as all steps within a given analysis
    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
        [ANALYSIS.RES.h_ttest(na, :), ANALYSIS.RES.holm_adjusted_alpha(na)] = multcomp_holm_bonferroni(ANALYSIS.RES.p_ttest(na,:), 'alpha', ANALYSIS.pstats); % Holm-Bonferroni correction      
        fprintf('The adjusted critical alpha for analysis %i is %1.6f   \n', na, ANALYSIS.RES.holm_adjusted_alpha(na));
    end % of for na loop

%__________________________________________________________________________    

case 3 % Strong FWER Control Permutation Test (Blaire-Karniski)
    
    fprintf('Performing corrections for multiple comparisons (permutation test)\n');
    
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
%__________________________________________________________________________    

case 4 % Cluster-Based Permutation Test
 
    fprintf('Performing corrections for multiple comparisons (cluster-based permutation test)\n');
    
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
%__________________________________________________________________________    

case 5 % KTMS Generalised FWER Control Using Permutation Testing
    
    fprintf('Performing corrections for multiple comparisons (KTMS generalised FWER control)\n');

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

%__________________________________________________________________________    

case 6 % Benjamini-Hochberg FDR Control

    fprintf('Performing corrections for multiple comparisons (Benjamini-Hochberg)\n');
    
    % Here a family of tests is defined as all steps within a given analysis
    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
        [ANALYSIS.RES.h_ttest(na, :), ANALYSIS.RES.bh_crit_alpha(na)] = multcomp_fdr_bh(ANALYSIS.RES.p_ttest(na,:), 'alpha', ANALYSIS.pstats);
        fprintf('The adjusted critical alpha for analysis %i is %1.6f \n', na, ANALYSIS.RES.bh_crit_alpha(na));
    end % of for na loop

%__________________________________________________________________________    

case 7 % Benjamini-Krieger-Yekutieli FDR Control
    
    fprintf('Performing corrections for multiple comparisons (Benjamini-Krieger-Yekutieli)\n');

    % Here a family of tests is defined as all steps within a given analysis
    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
        [ANALYSIS.RES.h_ttest(na, :), ANALYSIS.RES.bky_crit_alpha(na)] = multcomp_fdr_bky(ANALYSIS.RES.p_ttest(na,:), 'alpha', ANALYSIS.pstats);
        fprintf('The adjusted critical alpha for analysis %i is %1.6f \n', na, ANALYSIS.RES.bky_crit_alpha(na));
    end % of for na loop

%__________________________________________________________________________    

case 8 % Benjamini-Yekutieli FDR Control
    
    fprintf('Performing corrections for multiple comparisons (Benjamini-Yekutieli)\n');

    % Here a family of tests is defined as all steps within a given analysis
    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
        [ANALYSIS.RES.h_ttest(na, :), ANALYSIS.RES.by_crit_alpha(na)] = multcomp_fdr_by(ANALYSIS.RES.p_ttest(na,:), 'alpha', ANALYSIS.pstats);
        fprintf('The adjusted critical alpha for analysis %i is %1.6f \n', na, ANALYSIS.RES.by_crit_alpha(na));
    end % of for na loop

%__________________________________________________________________________    
% If some other option is chosen then do not correct, but notify user
otherwise
    fprintf('Unavailable multiple comparisons option chosen. Will use uncorrected p-values \n');
    ANALYSIS.RES.h_ttest = ANALYSIS.RES.h_ttest_uncorrected; 
end % of ANALYSIS.multcompstats switch