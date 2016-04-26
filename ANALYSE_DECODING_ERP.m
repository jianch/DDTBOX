function ANALYSE_DECODING_ERP(SLIST, sbjs_todo, dcg_todo, varargin)
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
% - SLIST (structure containing study configuration variables)
% - sbjs_todo (e.g., [1 2 3 4 6 7 9 10 13])
% - dcg_todo (discrimination group to analyse, as specified in SLIST.dcg_labels{dcg})
%
% optional:
% - allchan (analyse all possible channels? 1=yes / 2=no)
% - relchan (specify channels to be analysed (for temporal decoding only)
% - stmode (SPACETIME mode - 1=spatial / 2=temporal / 3=spatio-temporal)
% - avmode (AVERAGE mode - 1=no averaging; single-trial / 2=run average)
% - window_width_ms (width of sliding window in ms)
% - step_width_ms (step size with which sliding window is moved through the trial)
% - permstats (Testing against: 1=theoretical chance level / 2=permutation test results)
% - drawmode (Testing against: 1=average permutated distribution (default) / 2=random values drawn form permuted distribution (stricter))
% - pstats (critical p-value)
% - multcompstats (Correction for multiple comparisons:
% 0 = no correction
% 1 = Bonferroni correction
% 2 = Holm-Bonferroni correction
% 3 = Strong FWER Control Permutation Test
% 4 = Cluster-Based Permutation Test
% 5 = KTMS Generalised FWER Control
% 6 = Benjamini-Hochberg FDR Control
% 7 = Benjamini-Krieger-Yekutieli FDR Control
% 8 = Benjamini-Yekutieli FDR Control
% - nIterations (Number of permutation or bootstrap iterations for resampling-based multiple comparisons correction procedures)
% - KTMS_u (u parameter of the KTMS GFWER control procedure)
% - cluster_test_alpha (For cluster-based permutation test: Significance threshold for detecting effects at individual time windows)
% - disp_on (display a results figure? 0=no / 1=yes)
% - permdisp (display the results from permutation test in figure as separate line? 0=no / 1=yes)
% - disp_sign (display statistically significant steps in results figure? 0=no / 1=yes)
% - fw_do (analyse feature weights? 0=no / 1=yes)
% - fw_multcompstats (Feature weights correction for multiple comparisons:
% 1 = Bonferroni correction
% 2 = Holm-Bonferroni correction
% 3 = Strong FWER Control Permutation Test
% 4 = Cluster-Based Permutation Test (Currently not available)
% 5 = KTMS Generalised FWER Control
% 6 = Benjamini-Hochberg FDR Control
% 7 = Benjamini-Krieger-Yekutieli FDR Control
% 8 = Benjamini-Yekutieli FDR Control
%
% Options to display/not display feature weights results if they are analysed
% Averaged maps and stats:
% - fw_display_matrix (feature weights matrix. 0=no / 1=yes)
% - fw_display_average_zmap (z-standardised average feature weights. 0=no / 1=yes)
% - display_average_uncorr_threshmap (thresholded map uncorrected t-test results)
% - display_average_corr_threshmap (thresholded map of t-test results corrected for multiple comparisons)
% - fw_display_all_zmaps
% Individual maps and stats:
% - fw_display_all_zmaps (z-standardised average feature weights)
% - fw_display_all_uncorr_thresh_maps (thresholded map uncorrected t-test results)
% - fw_display_all_corr_thresh_maps (thresholded map of t-test results corrected for multiple comparisons)
%__________________________________________________________________________
%
% Variable naming convention: STRUCTURE_NAME.example_variable

%% GENERAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

% Define defaults at the beginning
options = struct(...
    'allchan', 1,...
    'relchan', [],...
    'stmode', 3,...
    'avmode', 1,...
    'window_width_ms', 10,...
    'step_width_ms', 10,...
    'permstats', 1,...
    'drawmode', 1,...
    'pstats', 0.05,...
    'multcompstats', 0,...
    'nIterations', 5000,...
    'KTMS_u', 1,...
    'cluster_test_alpha', 0.05,...
    'disp_on', 1,...
    'permdisp', 1,...
    'disp_sign', 1,...
    'fw_do', 1,...
    'fw_multcompstats', 1,...
    'fw_display_matrix', 1,...
    'fw_display_average_zmap', 0,...
    'fw_display_average_uncorr_threshmap', 0,...
    'fw_display_average_corr_threshmap', 0,...
    'fw_display_all_zmaps', 0,...
    'fw_display_all_uncorr_thresh_maps', 0,...
    'fw_display_all_corr_thresh_maps', 0 ...
    );

% Read the acceptable names
optionNames = fieldnames(options);

% Count arguments
nArgs = length(varargin);
if round(nArgs/2) ~= nArgs/2
   error([mfilename ' needs property name/property value pairs'])
end

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
   inpName = lower(pair{1}); % make case insensitive

   % Overwrite default options
   if any(strcmp(inpName,optionNames))
      options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name', inpName)
   end
end
clear pair
clear inpName

% Assignign variables to ANALYSIS structure
ANALYSIS.allchan = options.allchan; % Are all possible channels analysed? 1=yes (default if spatial/spatio-temporal) / 2=no
ANALYSIS.relchan = options.relchan; % specify channels to be analysed (for temporal only)

ANALYSIS.stmode = options.stmode; % SPACETIME mode (1=spatial / 2=temporal / 3=spatio-temporal)
ANALYSIS.avmode = options.avmode; % AVERAGE mode (1=no averaging; single-trial / 2=run average) 
ANALYSIS.window_width_ms = options.window_width_ms; % width of sliding window in ms
ANALYSIS.step_width_ms = options.step_width_ms; % step size with which sliding window is moved through the trial

ANALYSIS.permstats = options.permstats; % Testing against: 1=theoretical chance level / 2=permutation test results
ANALYSIS.drawmode = options.drawmode; % Testing against: 1=average permutated distribution (default) / 2=random values drawn form permuted distribution (stricter)

ANALYSIS.pstats = options.pstats; % critical p-value
ANALYSIS.multcompstats = options.multcompstats; % Correction for multiple comparisons: 
% 0 = no correction
% 1 = Bonferroni correction
% 2 = Holm-Bonferroni correction
% 3 = Strong FWER Control Permutation Test
% 4 = Cluster-Based Permutation Test
% 5 = KTMS Generalised FWER Control
% 6 = Benjamini-Hochberg FDR Control
% 7 = Benjamini-Krieger-Yekutieli FDR Control
% 8 = Benjamini-Yekutieli FDR Control
ANALYSIS.nIterations = options.nIterations; % Number of permutation or bootstrap iterations for resampling-based multiple comparisons correction procedures
ANALYSIS.KTMS_u = options.KTMS_u; % u parameter of the KTMS GFWER control procedure
ANALYSIS.cluster_test_alpha = options.cluster_test_alpha; % For cluster-based permutation test: Significance threshold for detecting effects at individual time windows (e.g. 0.05)

ANALYSIS.disp.on = options.disp_on; % display a results figure? 0=no / 1=yes
ANALYSIS.permdisp = options.permdisp; % display the results from permutation test in figure as separate line? 0=no / 1=yes
ANALYSIS.disp.sign = options.disp_sign; % display statistically significant steps in results figure? 0=no / 1=yes

ANALYSIS.fw.do = options.fw_do; % analyse feature weights? 0=no / 1=yes
ANALYSIS.fw.multcompstats = options.fw_multcompstats; % Feature weights correction for multiple comparisons:
% 1 = Bonferroni correction
% 2 = Holm-Bonferroni correction
% 3 = Strong FWER Control Permutation Test
% 4 = Cluster-Based Permutation Test (Currently not available)
% 5 = KTMS Generalised FWER Control
% 6 = Benjamini-Hochberg FDR Control
% 7 = Benjamini-Krieger-Yekutieli FDR Control
% 8 = Benjamini-Yekutieli FDR Control

% if feature weights are analysed, specify what is displayed
%__________________________________________________________________

% 0=no / 1=yes
ANALYSIS.fw.display_matrix = options.fw_display_matrix; % feature weights matrix

% averaged maps and stats
ANALYSIS.fw.display_average_zmap = options.fw_display_average_zmap; % z-standardised average FWs
ANALYSIS.fw.display_average_uncorr_threshmap = options.fw_display_average_uncorr_threshmap; % thresholded map uncorrected t-test results
ANALYSIS.fw.display_average_corr_threshmap = options.fw_display_average_corr_threshmap; % thresholded map of t-test results corrected for multiple comparisons

% individual maps and stats
ANALYSIS.fw.display_all_zmaps = options.fw_display_all_zmaps; % z-standardised average FWs
ANALYSIS.fw.display_all_uncorr_thresh_maps = options.fw_display_all_uncorr_thresh_maps; % thresholded map uncorrected t-test results
ANALYSIS.fw.display_all_corr_thresh_maps = options.fw_display_all_corr_thresh_maps; % thresholded map corrected t-test results (Bonferroni)

clear options;




% post-processing script = 2 - needed for interaction with other scripts to regulate
% functions such as saving data, calling specific sub-sets of parameters
global CALL_MODE
CALL_MODE = 3;

global DCGTODO;
DCGTODO = dcg_todo;



% define which subjects enter the second-level analysis
ANALYSIS.nsbj = size(sbjs_todo,2);
ANALYSIS.sbjs = sbjs_todo;
ANALYSIS.dcg_todo = dcg_todo;


fprintf('Group-level statistics will now be computed and displayed. \n'); 


%% OPEN FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

for s = 1:ANALYSIS.nsbj
    
    %% open subject data
    global SBJTODO;
    SBJTODO = s;
    sbj = ANALYSIS.sbjs(SBJTODO);
    
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
        
        % adjust for multiple comparisons 
        
        ANALYSIS.allsteps = size(ANALYSIS.data,2); % Calculate the number of tests (number of windows to analyse)
        
        % DF NOTE: The next 10 lines of code do not look like they are
        % related to multiple comparisons correction. Move to a different
        % section?
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
        
        if ANALYSIS.multcompstats == 1
            ANALYSIS.pstatsuse = ANALYSIS.pstats / ANALYSIS.allsteps; % Bonferroni correction
        elseif ANALYSIS.multcompstats ~= 1
            ANALYSIS.pstatsuse = ANALYSIS.pstats;
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
            [H,P] = ttest(ANALYSIS.RES.all_subj_acc(:,na,step),ANALYSIS.chancelevel,ANALYSIS.chancelevel); % simply against chance
            
        % test against permutation test results    
        elseif ANALYSIS.permstats == 2
            
            % against average permuted distribution
            if ANALYSIS.drawmode == 1
                
                [H,P] = ttest(ANALYSIS.RES.all_subj_acc(:,na,step),ANALYSIS.RES.all_subj_perm_acc(:,na,step),ANALYSIS.pstatsuse);
            
            % against one randomly drawn value (from all cross-val repetitions for each participant) for stricter test    
            elseif ANALYSIS.drawmode == 2
                
                for sbj = 1:ANALYSIS.nsbj
                    temp = randperm(size(ANALYSIS.RES.all_subj_perm_acc_reps_draw{sbj,na,step}(:,:),2));
                    drawone = temp(1); clear temp;
                    ANALYSIS.RES.draw_subj_perm_acc(sbj,na,step) = ANALYSIS.RES.all_subj_perm_acc_reps_draw{sbj,na,step}(1,drawone);
                    clear drawone;
                end % sbj
                
                [H,P] = ttest(ANALYSIS.RES.all_subj_acc(:,na,step),ANALYSIS.RES.draw_subj_perm_acc(:,na,step),ANALYSIS.pstatsuse);
                
            end % if ANALYSIS.drawmode
            
        end % if ANALYSIS.permstats
       
        ANALYSIS.RES.p_ttest(na,step) = P; clear P;
        ANALYSIS.RES.h_ttest(na,step) = H; clear H;
            
    end % of for step = 1:size(ANALYSIS.RES.mean_subj_acc,2) loop
    
end % of for na = 1:size(ANALYSIS.RES.mean_subj_acc,1) loop

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
