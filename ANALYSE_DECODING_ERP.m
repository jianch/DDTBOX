function ANALYSE_DECODING_ERP(study_name,vconf,input_mode,sbjs_todo,dcg_todo)
%__________________________________________________________________________
% DDTBOX script written by Stefan Bode 01/03/2013
%
% The toolbox was written with contributions from:
% Daniel Bennett, Jutta Stahl, Daniel Feuerriegel, Phillip Alday
%
% The author further acknowledges helpful conceptual input/work from: 
% Simon Lilburn, Philip L. Smith, Elaine Corbett, Carsten Murawski, 
% Carsten Bogler, John-Dylan Haynes
%__________________________________________________________________________
%
% This script is the master-script for the group-level analysis of EEG-decoding
% results. It will call several subscripts that run all possible analyses,
% depending on the specific decoding analyses.
%
% requires:
% - study_name (e.g. 'DEMO')
% - vconf (version of study configuration script, e.g., "1" for DEMO_config_v1)
% - input_mode (1 = use coded varialbles from first section / 2 = enter manually)
% - sbjs_todo (e.g., [1 2 3 4 6 7 9 10 13])
% - dcg_todo (discrimination group to analyse, as specified in SLIST.dcg_labels{dcg})

%__________________________________________________________________________
%
% Variable naming convention: STRUCTURE_NAME.example_variable

%% GENERAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

% post-processing script = 2 - needed for interaction with other scripts to regulate
% functions such as saving data, calling specific sub-sets of parameters
global CALL_MODE
CALL_MODE = 3;

global DCGTODO;
DCGTODO = dcg_todo;

sbj_list = [study_name '_config_v' num2str(vconf)]; % use latest slist-function!

% define which subjects enter the second-level analysis
ANALYSIS.nsbj = size(sbjs_todo,2);
ANALYSIS.sbjs = sbjs_todo;
ANALYSIS.dcg_todo = dcg_todo;

%% specify details about analysis & plotting

%__________________________________________________________________________
if input_mode == 0 % Hard-coded input

    % define all parameters of results to analyse & Plot
    %______________________________________________________________________
    ANALYSIS.allchan = 1; % Are all possible channels analysed? 1=yes (default if spatial/spatio-temporal) / 2=no
    ANALYSIS.relchan = []; % specify channels to be analysed (for temporal only)
        
    ANALYSIS.stmode = 3; % SPACETIME mode (1=spatial / 2=temporal / 3=spatio-temporal)
    ANALYSIS.avmode = 1; % AVERAGE mode (1=no averaging; single-trial / 2=run average) 
    ANALYSIS.window_width_ms = 10; % width of sliding window in ms
    ANALYSIS.step_width_ms = 10; % step size with which sliding window is moved through the trial
    
    ANALYSIS.permstats = 1; % Testing against: 1=theoretical chance level / 2=permutation test results
    ANALYSIS.drawmode = 1; % Testing against: 1=average permutated distribution (default) / 2=random values drawn form permuted distribution (stricter)
   
    ANALYSIS.pstats = 0.05; % critical p-value
    ANALYSIS.multcompstats = 4; % Correction for multiple comparisons: 
    % 0 = no correction
    % 1 = Bonferroni correction
    % 2 = Holm-Bonferroni correction
    % 3 = Strong FWER Control Permutation Test
    % 4 = Cluster-Based Permutation Test
    % 5 = KTMS Generalised FWER Control
    % 6 = Benjamini-Hochberg FDR Control
    % 7 = Benjamini-Krieger-Yekutieli FDR Control
    % 8 = Benjamini-Yekutieli FDR Control
    ANALYSIS.nIterations = 100; % Number of permutation or bootstrap iterations for resampling-based multiple comparisons correction procedures
    ANALYSIS.KTMS_u = 3; % u parameter of the KTMS GFWER control procedure
    ANALYSIS.cluster_test_alpha = 0.05; % For cluster-based test: Significance threshold for detecting effects at individual time windows (e.g. 0.05)
    
    ANALYSIS.disp.on = 1; % display a results figure? 0=no / 1=yes
    ANALYSIS.permdisp = 1; % display the results from permutation test in figure as separate line? 0=no / 1=yes
    ANALYSIS.disp.sign = 1; % display statistically significant steps in results figure? 0=no / 1=yes
    
    ANALYSIS.fw.do = 0; % analyse feature weights? 0=no / 1=yes
    
        % if feature weights are analysed, specify what is displayed
        %__________________________________________________________________
        
        % 0=no / 1=yes
        ANALYSIS.fw.display_matrix = 1; % feature weights matrix
        
        % averaged maps and stats
        ANALYSIS.fw.display_average_zmap = 0; % z-standardised average FWs
        ANALYSIS.fw.display_average_uncorr_threshmap = 0; % thresholded map uncorrected t-test results
        ANALYSIS.fw.display_average_corr_threshmap = 0; % thresholded map corrected t-test results (Bonferroni)
        
        % individual maps and stats
        ANALYSIS.fw.display_all_zmaps = 0; % z-standardised average FWs
        ANALYSIS.fw.display_all_uncorr_thresh_maps = 0; % thresholded map uncorrected t-test results
        ANALYSIS.fw.display_all_corr_thresh_maps = 0; % thresholded map corrected t-test results (Bonferroni)
%__________________________________________________________________________    

elseif input_mode == 1 % Prompted manual input
    
    % specify analysis channels
    ANALYSIS.allchan = input('Are all possible channels analysed? "0" for no; "1" for yes (default if spatial/spatio-temporal): ');
    
    if ANALYSIS.allchan ~= 1
        
        ANALYSIS.relchan = input('Enter the channels to be analysed (e.g. [1 4 5]): ');
        
    end
    
    % specify properties of the decoding analysis
    ANALYSIS.stmode = input('Specify the s/t-analysis mode of the original analysis. "1" spatial, "2" temporal. "3" spatio-temporal: ');
    ANALYSIS.avmode = input('Specify the average mode of the original analysis. "1" single-trial, "2" run-average: ');
    ANALYSIS.window_width_ms = input('Specify the window width [ms] of the original analysis: ');
    ANALYSIS.step_width_ms = input('Specify the step width [ms] of the original analysis: ');
    
    % specify stats
    ANALYSIS.permstats = input('Testing against: "1" chance level; "2" permutation distribution: ');
    
    if ANALYSIS.permstats == 2
        
        ANALYSIS.drawmode = input('Testing against: "1" average permutated distribution (default); "2" random values drawn form permuted distribution (stricter): ');
        ANALYSIS.permdisp = input('Do you wish to display chance-level test results in figure? "0" for no; "1" for yes: ');
        
    end
    
    ANALYSIS.pstats = input('Specify critical p-value for statistical testing (e.g. 0.05): ');
    ANALYSIS.multcompstats = input(['\nSpecify if you wish to control for multiple comparisons: \n"0" for no correction \n'...
        '"1" for Bonferroni \n"2" for Holm-Bonferroni \n"3" for Strong FWER Control Permutation Testing \n' ...
        '"4" for Cluster-Based Permutation Testing \n"5" for KTMS Generalised FWER Control \n' ...
        '"6" for Benjamini-Hochberg FDR Control \n"7" for Benjamini-Krieger-Yekutieli FDR Control \n' ...
        '"8" for Benjamini-Yekutieli FDR Control \n Option: ']);
    
    if ANALYSIS.multcompstats == 3 || ANALYSIS.multcompstats == 4 || ANALYSIS.multcompstats == 5
        ANALYSIS.nIterations = input('Number of bootstrap iterations for multiple comparisons procedure (at least 5000 is recommended): ');    
    end
    if ANALYSIS.multcompstats == 5
       ANALYSIS.KTMS_u = input('Enter the u parameter for the KTMS Generalised FWER control procedure: '); 
    end
    if ANALYSIS.multcompstats == 4
       ANALYSIS.cluster_test_alpha = input('Enter the significance threshold for detecting effects at individual time points (e.g. 0.05): '); 
    end
    
    % specify display options
    ANALYSIS.disp.on = input('Do you wish to display the results in figure(s)? "0" for no; "1" for yes: ');
    ANALYSIS.disp.sign = input('Specify if you wish to highlight significant results in figure. "0" for no; "1" for yes: ');
    
    % analyse feature weights
    ANALYSIS.fw.do = input('Do you wish to analyse the feature weights (only for spatial or spatio-temporal decoding)? "0" for no; "1" for yes: ');
    
    if ANALYSIS.fw.do == 1
        
        ANALYSIS.fw.display_average_zmap = input('Do you wish to display the group-level averaged, z-standardised feature weights as a heat map? "0" for no; "1" for yes: '); % z-standardised average FWs
        ANALYSIS.fw.display_average_uncorr_threshmap = input(...
            'Do you wish to display the statistical threshold map (uncorrected) for the group-level averaged, z-standardised feature weights as a heat map? "0" for no; "1" for yes: '); % thresholded map uncorrected t-test results
        ANALYSIS.fw.display_average_corr_threshmap = input(...
            'Do you wish to display the statistical threshold map (Bonferroni-corrected) for the group-level averaged, z-standardised feature weights as a heat map? "0" for no; "1" for yes: '); % thresholded map corrected t-test results (Bonferroni)
        
        % individual maps and stats
        ANALYSIS.fw.display_all_zmaps = input('');
        ANALYSIS.fw.display_all_uncorr_thresh_maps = input(...
            'Do you wish to display the statistical threshold map (uncorrected) for the group-level z-standardised feature weights for each time-step as a heat map? "0" for no; "1" for yes: ');
        ANALYSIS.fw.display_all_corr_thresh_maps = input(...
            'Do you wish to display the statistical threshold map (Bonferroni-corrected) for the group-level z-standardised feature weights for each time-step as a heat map? "0" for no; "1" for yes: ');
        
    end
    
end % input
%__________________________________________________________________________

fprintf('Group-level statistics will now be computed and displayed. \n'); 


%% OPEN FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

for s = 1:ANALYSIS.nsbj
    
    %% open subject data
    global SBJTODO;
    SBJTODO = s;
    sbj = ANALYSIS.sbjs(SBJTODO);
    
    global SLIST;
    eval(sbj_list);
    
    % open subject's decoding results       
    if size(dcg_todo,2) == 1
        
        fprintf('Loading results for subject %d in DCG %s.\n',sbj,SLIST.dcg_labels{dcg_todo});
        
        open_name = [(SLIST.output_dir) study_name '_SBJ' num2str(sbj) '_win' num2str(ANALYSIS.window_width_ms) '_steps' num2str(ANALYSIS.step_width_ms)...
            '_av' num2str(ANALYSIS.avmode) '_st' num2str(ANALYSIS.stmode) '_DCG' SLIST.dcg_labels{ANALYSIS.dcg_todo} '.mat'];

    elseif size(dcg_todo,2) == 2
        
        fprintf('Loading results for subject %d for cross decoding DCG %s => DCG %s.\n',sbj,SLIST.dcg_labels{dcg_todo(1)},SLIST.dcg_labels{dcg_todo(2)});
        
        open_name=[(SLIST.output_dir) study_name '_SBJ' num2str(sbj) '_win' num2str(ANALYSIS.window_width_ms) '_steps' num2str(ANALYSIS.step_width_ms)...
            '_av' num2str(ANALYSIS.avmode) '_st' num2str(ANALYSIS.stmode) '_DCG' SLIST.dcg_labels{ANALYSIS.dcg_todo(1)}...
            'toDCG' SLIST.dcg_labels{ANALYSIS.dcg_todo(2)} '.mat'];
    end   
   
    load(open_name);
    fprintf('Done.\n');
    
    ANALYSIS.analysis_mode=STUDY.analysis_mode;
    ANALYSIS.pointzero=SLIST.pointzero;
    
        
    %% fill in parameters and extract results 
    %______________________________________________________________________
    %
    % RESULTS contains averaged results:
    % RESULTS.subj_acc(analysis/channel,time-step) 
    % RESULTS.subj_perm_acc(analysis/channel,time-step) 
    % RESULTS contains raw results:
    % RESULTS.prediction_accuracy{analysis/channel}(time-step,cross-val_step,rep_step)
    % RESULTS.perm_prediction_accuracy{analysis/channel}(time-step,cross-val_step,rep_step)
    %
    % this section adds group results to ANALYSIS:
    % ANALYSIS.RES.all_subj_acc(subject,analysis/channel,time_step(fist_step:last_step))
    % ANALYSIS.RES.all_subj_perm_acc(subject,analysis/channel,time_step(fist_step:last_step))
    % ANALYSIS.RES.all_subj_perm_acc_reps(subject,analysis/channel,time_step(fist_step:last_step),cross-val_step,rep_step)
    
    % Define missing parameters using the first subject's dataset
    %______________________________________________________________________
    if s == 1 
        
        % ask for the specific time steps to analyse
        if ANALYSIS.avmode == 1 || ANALYSIS.avmode == 1 % DF NOTE: Is the second IF statement supposed to specify a different value?
    
            fprintf('\n');
            fprintf('You have %d time-steps in your RESULTS. Each time-step represents a %d ms time-window. \n',size(RESULTS.subj_acc,2),STUDY.window_width_ms);
            ANALYSIS.firststep = 1;
            ANALYSIS.laststep = input('Enter the number of the last time-window you want to analyse: ');

        end
    
        % shift everything back by step-width, as first bin gets label=0ms
        ANALYSIS.firststepms = (ANALYSIS.firststep * STUDY.step_width_ms) - STUDY.step_width_ms;
        ANALYSIS.laststepms = (ANALYSIS.laststep * STUDY.step_width_ms) - STUDY.step_width_ms;

        % create matrix for data indexing
        ANALYSIS.data(1,:) = 1:size(RESULTS.subj_acc,2); % for XTick
        ANALYSIS.data(2,:) = 0:STUDY.step_width_ms:( (size(RESULTS.subj_acc,2) - 1) * STUDY.step_width_ms); % for XLabel
        ptz = find(ANALYSIS.data(2,:) == ANALYSIS.pointzero); % find data with PointZero
        ANALYSIS.data(3,ptz) = 1; clear ptz; % for line location in plot

        % copy parameters from the config file
        ANALYSIS.step_width = STUDY.step_width;
        ANALYSIS.window_width = STUDY.window_width;
        ANALYSIS.sampling_rate = STUDY.sampling_rate;
        ANALYSIS.feat_weights_mode = STUDY.feat_weights_mode;
        
        ANALYSIS.nchannels = SLIST.nchannels;
                
        ANALYSIS.channellocs = SLIST.channellocs;
        ANALYSIS.channel_names_file = SLIST.channel_names_file;     
                
        % extract Tick/Labels for x-axis
        for datastep = 1:ANALYSIS.laststep
            ANALYSIS.xaxis_scale(1,datastep) = ANALYSIS.data(1,datastep);
            ANALYSIS.xaxis_scale(2,datastep) = ANALYSIS.data(2,datastep);
            ANALYSIS.xaxis_scale(3,datastep) = ANALYSIS.data(3,datastep);
        end
        
        % Define chance level for statistical analyses based on the
        % analysis type
        if STUDY.analysis_mode == 1 || STUDY.analysis_mode == 2
            ANALYSIS.chancelevel = ( 100 / size(SLIST.dcg{ANALYSIS.dcg_todo(1)},2) );
        elseif STUDY.analysis_mode == 3 || STUDY.analysis_mode == 4
            ANALYSIS.chancelevel = 0;
        end
        
        % Define channels to be used for group-analyses
        if ANALYSIS.allchan == 1

            % use all channels (default for spatial / spatial-temporal)
            ANALYSIS.allna = size(RESULTS.subj_acc,1);

        elseif ANALYSIS.allchan ~= 1

            % use specified number of channels
            ANALYSIS.allna = size(ANALYSIS.relchan,2);

        end
        
        
        if STUDY.analysis_mode == 1 || STUDY.analysis_mode == 2
            if size(ANALYSIS.dcg_todo,2) == 1
                ANALYSIS.DCG = SLIST.dcg_labels{ANALYSIS.dcg_todo};
            elseif size(ANALYSIS.dcg_todo,2) == 2
                ANALYSIS.DCG{1} = SLIST.dcg_labels{ANALYSIS.dcg_todo(1)};
                ANALYSIS.DCG{2} = SLIST.dcg_labels{ANALYSIS.dcg_todo(2)};
            end
        elseif STUDY.analysis_mode == 3 || STUDY.analysis_mode == 4    
            ANALYSIS.DCG = 'SVR_regression';
        end
        
                
    end % of if s == 1 statement
    
    %% extract results data from specified time-steps / channels
    %______________________________________________________________________
    
    for na = 1:ANALYSIS.allna
        
        % Extract classifier and permutation test accuracies
        ANALYSIS.RES.all_subj_acc(s,na,ANALYSIS.firststep:ANALYSIS.laststep) = RESULTS.subj_acc(na,ANALYSIS.firststep:ANALYSIS.laststep);
        ANALYSIS.RES.all_subj_perm_acc(s,na,ANALYSIS.firststep:ANALYSIS.laststep) = RESULTS.subj_perm_acc(na,ANALYSIS.firststep:ANALYSIS.laststep);
            
        % needed if one wants to test against distribution of randomly
        % drawn permutation results (higher variance, stricter testing)
        ANALYSIS.RES.all_subj_perm_acc_reps(s,na,ANALYSIS.firststep:ANALYSIS.laststep,:,:) = RESULTS.perm_prediction_accuracy{na}(ANALYSIS.firststep:ANALYSIS.laststep,:,:);
            
    end
    %______________________________________________________________________
    
    % Extract feature weights
    if ~isempty(RESULTS.feature_weights)
        ANALYSIS.RES.feature_weights{s} = RESULTS.feature_weights{1};
    end
    
    clear RESULTS;
    clear STUDY;
    
end % of for n = 1:ANALYSIS.nsbj loop

fprintf('All data from all subjects loaded.\n');

%% AVERAGE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

% Calculate average accuracy & standard error across subjects
M(:,:) = mean(ANALYSIS.RES.all_subj_acc,1);
ANALYSIS.RES.mean_subj_acc(:,:) = M'; clear M;

SE(:,:) = (std(ANALYSIS.RES.all_subj_acc,1))/(sqrt(ANALYSIS.nsbj));
ANALYSIS.RES.se_subj_acc(:,:) = SE'; clear SE;

if ANALYSIS.permstats == 2
    
    % OPTION 1: Use average results from random-labels test
    % Calculate average accuracy & standard error across subjects for permutation results
    M(:,:) = mean(ANALYSIS.RES.all_subj_perm_acc,1);
    ANALYSIS.RES.mean_subj_perm_acc(:,:) = M'; clear M;
    
    SE(:,:) = (std(ANALYSIS.RES.all_subj_perm_acc,1)) / (sqrt(ANALYSIS.nsbj));
    ANALYSIS.RES.se_subj_perm_acc(:,:) = SE'; clear SE;

    % OPTION 2: draw values from random-labels test
    % average permutation results across cross-validation steps, but draw later 
    % one for each participant for statistical testing!
    for subj = 1:ANALYSIS.nsbj
        for ana = 1:ANALYSIS.allna
            for step = 1:ANALYSIS.laststep
                temp(:,:) = ANALYSIS.RES.all_subj_perm_acc_reps(subj,ana,step,:,:);
                mtemp = mean(temp,1);
                ANALYSIS.RES.all_subj_perm_acc_reps_draw{subj,ana,step} = mtemp;
                clear temp; clear mtemp;
            end % step
        end % ana
    end % sbj

end % of if ANALYSIS.permstats == 2 statement

fprintf('All data from all subjects averaged.\n');

%% STATISTICAL TESTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

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

% Checking if an inappropriate multiple comparisons correction procedure
% has been selected for paired-samples comparisons against permutation
% distribution.
if ANALYSIS.multcompstats == 3 || ANALYSIS.multcompstats == 4 || ANALYSIS.multcompstats == 5
    if ANALYSIS.permstats == 2; % If user chooses to compare against permutation distribution instead of against chance level.
        fprintf('This multiple comparisons correction procedure cannot be used for comparison against permutation results. Switching to Bonferroni correction...\n');
        ANALYSIS.multcompstats = 1;
    end     
end

fprintf('Performing corrections for multiple comparisons...\n');

ANALYSIS.RES.h_ttest = zeros(size(ANALYSIS.RES.mean_subj_acc,1), size(ANALYSIS.RES.mean_subj_acc,2));
n_total_steps = size(ANALYSIS.RES.mean_subj_acc,2); % Moving to this variable for easier interpretation of the code

switch ANALYSIS.multcompstats

case 0 % No correction for multiple comparisons

    ANALYSIS.RES.h_ttest = ANALYSIS.RES.h_ttest_uncorrected; 

%__________________________________________________________________________

case 1 % Bonferroni Correction

ANALYSIS.pstats_bonferroni_corrected = ANALYSIS.pstats / n_total_steps; % Bonferroni correction

for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
    for step = 1:n_total_steps % step
        % Test against Bonferroni-corrected threshold
        if ANALYSIS.RES.p_ttest(na,step) < ANALYSIS.pstats_bonferroni_corrected
            ANALYSIS.RES.h_ttest(na,step) = 1;
        end
    end
end

%__________________________________________________________________________    

case 2 % Holm-Bonferroni Correction

% Here a family of tests is defined as all steps within a given analysis
for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
    % Sort p-values from smallest to largest
    sorted_p = sort(ANALYSIS.RES.p_ttest(na,:));
    foundCritAlpha = 0; % Reset to signify we have not found Holm-corrected critical alpha level

    % Find critical alpha level (sig. for all smaller p-values than this)
    for holmStep = 1:n_total_steps
       if sorted_p(holmStep) > ANALYSIS.pstats / (n_total_steps + 1 - holmStep) & foundCritAlpha == 0
           ANALYSIS.holm_corrected_alpha(na) = sorted_p(holmStep);
           foundCritAlpha = 1;
       end
    end

    if ~exist('ANALYSIS.holm_corrected_alpha', 'var') % If all null hypotheses are rejected
        ANALYSIS.holm_corrected_alpha(na) = 0;
    end

    % Declare tests significant if they are smaller than the adjusted critical alpha
    for step = 1:n_total_steps
        if ANALYSIS.RES.p_ttest(na, step) < ANALYSIS.holm_corrected_alpha(na)   
            ANALYSIS.RES.h_ttest(na,step) = 1;
        else
            ANALYSIS.RES.h_ttest(na,step) = 0;
        end
    end      
end % of for na loop

%__________________________________________________________________________    

case 3 % Strong FWER Control Permutation Test

    % Note: this correction will only work if the user has run
    % permutation tests during the decoding stage.

    % Check if permutation decoding results have been obtained, and if not,
    % warn the user and skip MC correction.
    
    % TODO: write this part. Process an ID without doing permutation
    % testing and see what the results would look like.
    
    
    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
    
    % Generate t(max) distribution from the null data (use permutation results
    % when implementing in DDTBOX)
    ANALYSIS.t_max(na, 1:ANALYSIS.nIterations) = zeros(1, ANALYSIS.nIterations);
    t_stat(1:n_total_steps, 1:ANALYSIS.nIterations) = zeros(n_total_steps, ANALYSIS.nIterations);
        for iteration = 1:ANALYSIS.nIterations
            clear temp; % Clearing out temp variable
            % Draw a random sample for each test
            for step = 1:n_total_steps
                % Draw a bootstrap sample (i.e. sampling with replacement)
                temp(1:ANALYSIS.nsbj, step) = randsample(ANALYSIS.RES.all_subj_perm_acc(1:ANALYSIS.nsbj,na, step),  ANALYSIS.nsbj, true);  
                [~, ~, ~, temp_stats] = ttest(temp(1:ANALYSIS.nsbj, step), 0, 'Alpha', ANALYSIS.pstats);
                t_stat(step, iteration) = abs(temp_stats.tstat);
            end    

            % Get the maximum t-value within the family of tests and store in a
            % vector. This is to create a null hypothesis distribution.
            ANALYSIS.t_max(na, iteration) = max(t_stat(:, iteration));  
        end % of for iteration loop
    
        % Calculating the 95th percentile of t_max values (used as decision
        % critieria for statistical significance)
        ANALYSIS.permtest_null_cutoff(na) = prctile(ANALYSIS.t_max(na, 1:ANALYSIS.nIterations), ((1 - ANALYSIS.pstats) * 100));

        % Checking whether each test statistic is above the specified threshold:
        for step = 1:n_total_steps
            if abs(ANALYSIS.RES.t_ttest(na,step)) > ANALYSIS.permtest_null_cutoff(na)
                ANALYSIS.RES.h_ttest(na,step) = 1;
            else
                ANALYSIS.RES.h_ttest(na,step) = 0;
            end
        end % of for step loop
        
    end % of for na loop
    
%__________________________________________________________________________    

case 4 % Cluster-Based Permutation Test
 
    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
    
            % Generate a maximum cluster mass distribution from the permutation
            % test results
            max_cluster_mass = zeros(1, ANALYSIS.nIterations);
            cluster_perm_test_h = zeros(n_total_steps, ANALYSIS.nIterations);
            clear t_stat;

            for iteration = 1:ANALYSIS.nIterations;
                % Draw a random bootstrap sample for each test
                for step = 1:n_total_steps  
                    % Draw a bootstrap sample (i.e. sampling with replacement)
                    temp(1:ANALYSIS.nsbj, step) = randsample(ANALYSIS.RES.all_subj_perm_acc(1:ANALYSIS.nsbj,na, step),  ANALYSIS.nsbj, true);

                    [cluster_perm_test_h(step, iteration), ~, ~, temp_stats] = ttest(temp(1:ANALYSIS.nsbj, step), 0, 'Alpha', ANALYSIS.cluster_test_alpha);
                    t_stat(step, iteration) = temp_stats.tstat; % Get t statistic
                    % Marking the sign of each t statistic to avoid clustering pos
                    % and neg significant results
                    if t_stat(step, iteration) < 0;
                        t_sign(step, iteration) = -1; 
                    else
                        t_sign(step, iteration) = 1; 
                    end
                end    

                % Identify clusters and generate a maximum cluster statistic
                cluster_mass_vector = [0]; % Resets vector of cluster masses
                cluster_counter = 0;

                for step = 1:n_total_steps     
                    if cluster_perm_test_h(step, iteration) == 1
                        if step == 1 % If the first test in the set
                            cluster_counter = cluster_counter + 1;
                            cluster_mass_vector(cluster_counter) = abs(t_stat(step, iteration));
                        else
                            % Add to the cluster if there are consecutive
                            % statistically significant tests with the same sign.
                            % Otherwise, make a new cluster.
                            if cluster_perm_test_h(step - 1, iteration) == 1 && t_sign(step - 1, iteration) == t_sign(step, iteration)
                                cluster_mass_vector(cluster_counter) = cluster_mass_vector(cluster_counter) + abs(t_stat(step, iteration));
                            else
                                cluster_counter = cluster_counter + 1;
                                cluster_mass_vector(cluster_counter) = abs(t_stat(test, iteration));
                            end 
                        end % of if test == 1
                    end % of if clusterPermTest
                end % of for steps = 1:n_total_steps

                % Find the maximum cluster mass
                max_cluster_mass(iteration) = max(cluster_mass_vector);
            end % of iterations loop

            % Calculating the 95th percentile of maximum cluster mass values (used as decision
            % critieria for statistical significance)
            cluster_mass_null_cutoff = prctile(max_cluster_mass(iteration), ((1 - ANALYSIS.pstats) * 100));


            % Calculate cluster masses in the actual (non-permutation) tests
            cluster_mass_vector = [0]; % Resets vector of cluster masses
            cluster_counter = 0;
            cluster_locations = zeros(1, n_total_steps);
            cluster_corrected_sig_steps = zeros(1, n_total_steps);
            clear t_sign;

            for step = 1:n_total_steps   
                if ANALYSIS.RES.h_ttest_uncorrected(na,step) == 1
                    if step == 1 % If the first test in the set
                        cluster_counter = cluster_counter + 1;
                        cluster_mass_vector(cluster_counter) = abs(ANALYSIS.RES.t_ttest(na,step));
                        cluster_locations(step) = cluster_counter;
                        % Tagging as positive or negative sign effect
                        if ANALYSIS.RES.t_ttest(na,step) < 0
                            t_sign(step) = -1;
                        else
                            t_sign(step) = 1;
                        end
                elseif step > 1
                    % Tagging as positive or negative sign effect
                    if ANALYSIS.RES.t_ttest(na,step) < 0
                        t_sign(step) = -1;
                    else
                        t_sign(step) = 1;
                    end

                    % Add to the same cluster only if the previous test was sig.
                    % and of the same sign (direction).
                    if ANALYSIS.RES.h_ttest_uncorrected(na,step - 1) == 1 && t_sign(step - 1) == t_sign(step)
                        cluster_mass_vector(cluster_counter) = cluster_mass_vector(cluster_counter) + abs(ANALYSIS.RES.t_ttest(na,step));
                        cluster_locations(step) = cluster_counter;
                    else
                        cluster_counter = cluster_counter + 1;
                        cluster_mass_vector(cluster_counter) = abs(ANALYSIS.RES.t_ttest(na,step));
                        cluster_locations(step) = cluster_counter;
                    end 
                end % of if step == 1
            end % of if ANALYSIS.RES.h_ttest_uncorrected(na,step) == 1  
        end % of for step = 1:n_total_steps

        for cluster_no = 1:length(cluster_mass_vector);
            if cluster_mass_vector(cluster_no) > cluster_mass_null_cutoff
                cluster_corrected_sig_steps(cluster_locations == cluster_no) = 1;
            end
        end

        % Update analysis structure with cluster-corrected significant time
        % windows
        ANALYSIS.RES.h_ttest(na, :) = cluster_corrected_sig_steps;
      
        
    end % of na loop  

%__________________________________________________________________________    

case 5 % KTMS Generalised FWER Control Using Permutation Testing
    
   
    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
        
        % Make a vector to denote statistically significant steps
        KTMS_sig_effect_locations = zeros(1, n_total_steps);

        % Sort p-values from smallest to largest
        sorted_p = sort(ANALYSIS.RES.p_ttest(na,:));
        
        % Automatically reject the u smallest hypotheses (u is set by user as KTMS_u variable).
        KTMS_autoRejectAlpha = sorted_p(ANALYSIS.KTMS_u);
        KTMS_sig_effect_locations(ANALYSIS.RES.p_ttest(na, :) <= KTMS_autoRejectAlpha) = 1; % Mark tests with u smallest p-values as statistically significant.
        
        % Run strong FWER control permutation test but use u + 1th most extreme
        % test statistic.
        KTMS_t_max = zeros(1, ANALYSIS.nIterations);
        t_stat = zeros(n_total_steps, ANALYSIS.nIterations);
        clear temp; % Clear temp variable
        
        for iteration = 1:ANALYSIS.nIterations

            % Draw a random sample for each test
            for step = 1:n_total_steps
                % Draw a bootstrap sample (i.e. sampling with replacement)
                temp(1:ANALYSIS.nsbj, step) = randsample(ANALYSIS.RES.all_subj_perm_acc(1:ANALYSIS.nsbj,na, step),  ANALYSIS.nsbj, true);  
                [~, ~, ~, temp_stats] = ttest(temp(1:ANALYSIS.nsbj, step), 0, 'Alpha', ANALYSIS.pstats);
                t_stat(step, iteration) = abs(temp_stats.tstat);
            end    

            % Get the maximum t-value within the family of tests and store in a
            % vector. This is to create a null hypothesis distribution.
            t_sorted = sort(t_stat(:, iteration), 'descend');
            KTMS_t_max(iteration) = t_sorted(ANALYSIS.KTMS_u + 1);
        end

        % Calculating the 95th percentile of t_max values (used as decision
        % critieria for statistical significance)
        ANALYSIS.KTMS_Null_Cutoff(na) = prctile(KTMS_t_max, ((1 - ANALYSIS.pstats) * 100));

        % Checking whether each test statistic is above the specified threshold:
        for step = 1:n_total_steps
            if abs(ANALYSIS.RES.t_ttest(na,step)) > ANALYSIS.KTMS_Null_Cutoff(na);
                KTMS_sig_effect_locations(step) = 1;
            end
        end
        
        % Marking statistically significant tests in the ANALYSIS structure
        ANALYSIS.RES.h_ttest(na,:) = KTMS_sig_effect_locations;    
    end    
    
%__________________________________________________________________________    

case 6 % Benjamini-Hochberg FDR Control

% Here a family of tests is defined as all steps within a given analysis
for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis
    % Sort p-values from smallest to largest
    sorted_p = sort(ANALYSIS.RES.p_ttest(na,:));

    % Find critical k value
    for benhoch_step = 1:n_total_steps
        if sorted_p(benhoch_step) <= (benhoch_step / n_total_steps) * ANALYSIS.pstats
            ANALYSIS.benhoch_critical_alpha(na) = sorted_p(benhoch_step);
        end
    end
    % If no steps are significant set critical alpha to zero
    if ~exist('benhoch_critical_alpha', 'var')
        ANALYSIS.benhoch_critical_alpha(na) = 0;
    end

    % Declare tests significant if they are smaller than or equal to the adjusted critical alpha
    for step = 1:n_total_steps
        if ANALYSIS.RES.p_ttest(na, step) <= ANALYSIS.benhoch_critical_alpha(na)   
            ANALYSIS.RES.h_ttest(na,step) = 1;
        else
            ANALYSIS.RES.h_ttest(na,step) = 0;
        end
    end     
end % of for na loop

%__________________________________________________________________________    

case 7 % Benjamini-Krieger-Yekutieli FDR Control

    % Here a family of tests is defined as all steps within a given analysis
    for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis

        % Stage 1: Estimate the number of false null hypotheses using the modified
        % alpha level.

        % Sort p-values from smallest to largest
        sorted_p = sort(ANALYSIS.RES.p_ttest(na,:));

        % Find critical k value (see tutorial notes)
        for BKY_step = 1:n_total_steps
            if sorted_p(BKY_step) <= (BKY_step / n_total_steps) * (ANALYSIS.pstats / ( 1 + ANALYSIS.pstats));
                BKY_stage1_critical_alpha = sorted_p(BKY_step);
            end
        end
        % If no tests are significant set critical alpha to zero
        if ~exist('BKY_stage1_critical_alpha', 'var')
            BKY_stage1_critical_alpha = 0;
        end

        % Declare tests significant if they are smaller than or equal to the adjusted critical alpha
        BKY_stage1_h = zeros(1, n_total_steps); % Preallocate for speed
        for step = 1:n_total_steps
            if ANALYSIS.RES.p_ttest(na, step) <= BKY_stage1_critical_alpha   
                BKY_stage1_h(step) = 1;
            else
                BKY_stage1_h(step) = 0;
            end
        end

        % Count the number of rejected null hypotheses (for use in step 2)
        BKY_stage1_n_rejections = sum(BKY_stage1_h);

        if BKY_stage1_n_rejections == 0; % if no null hypotheses were rejected
            ANALYSIS.RES.h_ttest(na,1:n_total_steps) = 0; % Don't reject any hypotheses

        elseif BKY_stage1_n_rejections == n_total_steps; % if all null hypotheseses were rejected
            ANALYSIS.RES.h_ttest(na,1:n_total_steps) = 1; % Reject all hypotheses

        else % If some (but not all) null hypotheses were rejected  
            for step = 1:n_total_steps
                if sorted_p(step) <= (step / n_total_steps) * ( (n_total_steps / (n_total_steps - BKY_stage1_n_rejections) ) * (ANALYSIS.pstats / ( 1 + ANALYSIS.pstats)) );
                    BKY_stage2_critical_alpha = sorted_p(step);
                end
            end

            % If no tests are significant set critical alpha to zero
            if ~exist('BKY_stage2_critical_alpha', 'var')
                BKY_stage2_critical_alpha = 0;
            end

            % Declare tests significant if they are smaller than or equal to the adjusted critical alpha
            for step = 1:n_total_steps
                if ANALYSIS.RES.p_ttest(na, step) <= BKY_stage2_critical_alpha   
                    ANALYSIS.RES.h_ttest(na,1:step) = 1;
                else
                    ANALYSIS.RES.h_ttest(na,1:step) = 0;
                end
            end

        end % of if BKY_stage1_n_rejections

    end % of for na loop

%__________________________________________________________________________    

case 8 % Benjamini-Yekutieli FDR Control

% Here a family of tests is defined as all steps within a given analysis
for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) % analysis

    % Sort p-values from smallest to largest
    sorted_p = sort(ANALYSIS.RES.p_ttest(na,:));

    % j values precalculated to help calculate the Benjamini-Yekutieli critical alpha
    j_values = zeros(1,n_total_steps);
    for j_iteration = 1:n_total_steps
        j_values(j_iteration) = 1 / j_iteration;
    end

    % Find critical k value (see tutorial notes)
    for benyek_step = 1:n_total_steps

        if sorted_p(benyek_step) <= (benyek_step / n_total_steps * sum(j_values)) * ANALYSIS.pstats
            benyek_critical_alpha = sorted_p(benyek_step);
        end
    end
    % If no tests are significant set critical alpha to zero
    if ~exist('benyek_critical_alpha', 'var')
        benyek_critical_alpha = 0;
    end

    % Declare tests significant if they are smaller than or equal to the adjusted critical alpha
    for step = 1:n_total_steps
        if ANALYSIS.RES.p_ttest(na, step) <= benyek_critical_alpha             
            ANALYSIS.RES.h_ttest(na,step) = 1;
        else
            ANALYSIS.RES.h_ttest(na,step) = 0;
        end
    end
end % of for na loop

%__________________________________________________________________________    
% If some other option is chosen then do not correct, but notify user
otherwise
    fprintf('Unavailable multiple comparisons option chosen. Will use uncorrected p-values \n');
    ANALYSIS.RES.h_ttest = ANALYSIS.RES.h_ttest_uncorrected; 
end % of ANALYSIS.multcompstats switch


fprintf('All group statistics performed.\n');

%% FEATURE WEIGHT ANALYSIS
%__________________________________________________________________________

if ANALYSIS.fw.do == 1
    
    [FW_ANALYSIS] = analyse_feature_weights_erp(ANALYSIS);
    
else
    
    FW_ANALYSIS = [];
    
end


%% SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

if size(dcg_todo,2) == 1 % Standard decoding analyses

    savename = [(SLIST.output_dir) study_name '_GROUPRES_NSBJ' num2str(ANALYSIS.nsbj) '_win' num2str(ANALYSIS.window_width_ms) '_steps' num2str(ANALYSIS.step_width_ms)...
        '_av' num2str(ANALYSIS.avmode) '_st' num2str(ANALYSIS.stmode) '_DCG' SLIST.dcg_labels{ANALYSIS.dcg_todo} '.mat'];
    
elseif size(dcg_todo,2) == 2 % Cross-condition decoding analyses
    
    savename = [(SLIST.output_dir) study_name '_GROUPRES_NSBJ' num2str(ANALYSIS.nsbj) '_win' num2str(ANALYSIS.window_width_ms) '_steps' num2str(ANALYSIS.step_width_ms)...
        '_av' num2str(ANALYSIS.avmode) '_st' num2str(ANALYSIS.stmode) '_DCG' SLIST.dcg_labels{ANALYSIS.dcg_todo(1)}...
        'toDCG' SLIST.dcg_labels{ANALYSIS.dcg_todo(2)} '.mat'];

end

save(savename,'ANALYSIS','FW_ANALYSIS');

fprintf('All results saved in %s. \n',savename);


%% PLOT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

if ANALYSIS.disp.on == 1
    
    fprintf('Results will be plotted. \n');
    display_group_results_erp(ANALYSIS);
    
elseif ANALYSIS.disp.on ~= 1
    
    fprintf('No figures were produced for the results. \n');
    
end

%__________________________________________________________________________
