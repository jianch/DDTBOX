% This script contains the SLIST config settings and also each variadic
% input that can be input into DECODING_ERP and ANALYSE_DECODING_ERP. Each
% set of variadic inputs (and an example call for each function) are
% inserted below.



%% SLIST configuration settings
% This is input as the SLIST structure to DECODING_ERP and ANALYSE_DECODING_ERP
%
% All copied from the config file DEMO_config_v1

bdir = 'F:\MVPA_WORKSHOP\'; % Base directory
output_dir = 'F:\MVPA_WORKSHOP\DECODING_RESULTS\preanalysed\'; % Directory in which the decoding results will be saved


% Filepaths and filenames of the EEG data files
sbj_code = {...

    ['DATA\sbj1\SBJ1_full'];... %1
    ['DATA\sbj2\SBJ2_full'];... %2 
    ['DATA\sbj3\SBJ3_full'];... %3
    ['DATA\sbj4\SBJ4_full'];... %4
    ['DATA\sbj5\SBJ5_full'];... %5

    };

% subject parameters
SLIST.number = sn;
SLIST.sbj_code = sbj_code{sn};
SLIST.output_dir = output_dir;
SLIST.data_struct_name = 'eeg_sorted_cond';

% channels    
SLIST.nchannels = 64; % Number of channels in the dataset
SLIST.channels = 'channel_labels'; 
SLIST.channel_names_file = 'channel_inf.mat'; % Name of the .mat file containing channel information
SLIST.channellocs = [bdir 'locations\']; % Directory of the .mat file containing channel information
SLIST.eyes = []; % Channel indices of ocular electrodes
SLIST.extra = [0]; % Channel indices of electrodes to exclude from the classification analyses

% sampling rate and baseline
SLIST.sampling_rate = 1000; % Sampling rate (Hz)
SLIST.pointzero = 100; % Corresponds to time zero, for example stimulus onset (in ms, from the beginning of the epoch)

% CREATE DCGs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

% Label each condition
% Example: SLIST.cond_labels{condition number} = 'Name of condition';
SLIST.cond_labels{1} = '0same';
SLIST.cond_labels{2} = '180same';
SLIST.cond_labels{3} = '0diff';
SLIST.cond_labels{4} = '180diff';

% Discrimination groups
% Enter the condition numbers of the conditions to discriminate between
% Example: SLIST.dcg{Discrimination group number} = [condition number 1, condition number 2];
SLIST.dcg{1} = [1 3]; % 
SLIST.dcg{2} = [2 4]; % 

% Label each discrimination group
% Example: SLIST.dcg_labels{Discrimination group number} = 'Name of discrimination group'
SLIST.dcg_labels{1} = '0SD';
SLIST.dcg_labels{2} = '180SD';

SLIST.ndcg = size(SLIST.dcg,2);
SLIST.nclasses = size(SLIST.dcg{1},2);      
SLIST.ncond = size(SLIST.cond_labels,2);

SLIST.data_open_name = [bdir (sbj_code{sn}) '.mat'];
SLIST.data_save_name = [bdir (sbj_code{sn}) '_data.mat'];








%% DECODING_ERP Variadic inputs

if input_mode == 0 % Hard-coded input

    % get study parameters (hardcoded)
    stmode = 3; % SPACETIME mode (1=spatial / 2=temporal / 3=spatio-temporal)
    avmode = 1; % AVERAGE mode (1=no averaging; single-trial / 2=run average) 
    analysis_mode = 1;% ANALYSIS mode (1=SVM classification, 2=LDA classification, 3=SVR increments, 4=SVR continous)
    
    window_width_ms = 10; % width of sliding window in ms
    step_width_ms = 10; % step size with which sliding window is moved through the trial
    
    perm_test = 1; % run the permutation-decoding? 0=no / 1=yes
    perm_disp = 1; % display the permutation results in figure? 0=no / 1=yes
    display_on = 1; % 1=figure displayed, 0=no figure
    
    rt_match = 0; % Use RT-matching algorithm to select trials? 0=no / 1=yes
    feat_weights_mode = 1; % Extract feature weights? 0=no / 1=yes
    cross_val_steps = 10; % How many cross-validation steps (if no runs available)?
    n_rep_cross_val = 10; % How many repetitions of full cross-validation with randomly re-ordered data?
    permut_rep = 10; % How many repetitions of full cross-validation with permutation results?

elseif input_mode == 1 % Prompt user for input

    % get study parameters via input
    stmode = input('Enter "1" for spatial, "2" for temporal, or "3" for spatio-temporal decoding: ');
    avmode = input('Enter "1" for single-trial decoding or "2" for run-average decoding: ');
    if avmode == 1 % Options for single-trial decoding
        cross_val_steps = input('How many cross-validation steps do you wish to perform?');
        n_rep_cross_val = input('How many independent repetitions of the analysis do you wish to perform?');
        rt_match = input('Do you wish to RT-match the trials from the conditions to be decoded? (1=yes, 0=no) ');
    end
    analysis_mode = input('Specifiy analysis method: "1" for Class SVM, "2" for Class LDA, "3" increments SVR, "4" continuous SVR: '); 
    window_width_ms = input('Enter decoding window width in ms: ');
    step_width_ms = input('Enter step width for moving the decoding window in ms: ');
    cross = input('Do you wish to perform cross-condition decoding? "0" for no, "1" for yes:');
    if cross > 0
        dcgs = input('Enter two discriminations groups for cross-decoding (e.g.[1 2]):');
        dcg_todo = dcgs;
    end
        perm_test = input('Do you want to run a permutation test? (1=yes, 0=no) '); 
    if perm_test == 1
        permut_rep = input('How many repetitions of the permutation test do you wish to perform? ');
        perm_disp = input('Do you wish to plot the permutation test results? (1=yes, 0=no) '); 
    end
    display_on = input('Do you wish to plot the individual decoding results? (1=yes, 0=no) ');

end

% Subject, dcg_todo inputs:
sbj = 1; % Subject number to decode
dcg_todo = 1; % Discrimination group to decode
cross = 0; % Cross-condition decoding? 1=yes / 0=no

% Call to function with all variadic inputs:
DECODING_ERP(SLIST, sbj, dcg_todo, ...
    'stmode', stmode,...
    'avmode', avmode,...
    'analysis_mode', analysis_mode,...
    'window_width_ms', window_width_ms,...
    'step_width_ms', step_width_ms,...
    'cross', cross,...
    'perm_test', perm_test,...
    'perm_disp', perm_disp,...
    'display_on', display_on,...
    'rt_match', rt_match,...
    'feat_weights_mode', feat_weights_mode,...
    'cross_val_steps', cross_val_steps,...
    'n_rep_cross_val', n_rep_cross_val,...
    'permut_rep', permut_rep);











%% ANALYSE_DECODING_ERP variadic inputs

if input_mode == 0 % Hard-coded input

    % define all parameters of results to analyse & Plot
    %______________________________________________________________________
    allchan = 1; % Are all possible channels analysed? 1=yes (default if spatial/spatio-temporal) / 2=no
    relchan = []; % specify channels to be analysed (for temporal only)
        
    stmode = 3; % SPACETIME mode (1=spatial / 2=temporal / 3=spatio-temporal)
    avmode = 1; % AVERAGE mode (1=no averaging; single-trial / 2=run average) 
    window_width_ms = 10; % width of sliding window in ms
    step_width_ms = 10; % step size with which sliding window is moved through the trial
    
    permstats = 1; % Testing against: 1=theoretical chance level / 2=permutation test results
    drawmode = 1; % Testing against: 1=average permutated distribution (default) / 2=random values drawn form permuted distribution (stricter)
   
    pstats = 0.05; % critical p-value
    multcompstats = 0; % Correction for multiple comparisons: 
    % 0 = no correction
    % 1 = Bonferroni correction
    % 2 = Holm-Bonferroni correction
    % 3 = Strong FWER Control Permutation Test
    % 4 = Cluster-Based Permutation Test
    % 5 = KTMS Generalised FWER Control
    % 6 = Benjamini-Hochberg FDR Control
    % 7 = Benjamini-Krieger-Yekutieli FDR Control
    % 8 = Benjamini-Yekutieli FDR Control
    nIterations = 5000; % Number of permutation or bootstrap iterations for resampling-based multiple comparisons correction procedures
    KTMS_u = 2; % u parameter of the KTMS GFWER control procedure
    cluster_test_alpha = 0.05; % For cluster-based test: Significance threshold for detecting effects at individual time windows (e.g. 0.05)
    
    disp.on = 1; % display a results figure? 0=no / 1=yes
    permdisp = 1; % display the results from permutation test in figure as separate line? 0=no / 1=yes
    disp.sign = 1; % display statistically significant steps in results figure? 0=no / 1=yes
    
    fw_do = 0; % analyse feature weights? 0=no / 1=yes
    fw_multcompstats = 1; % Feature weights correction for multiple comparisons:
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
        fw_display_matrix = 1; % feature weights matrix % 0=no / 1=yes
        
        % averaged maps and stats
        fw_display_average_zmap = 0; % z-standardised average FWs
        fw_display_average_uncorr_threshmap = 0; % thresholded map uncorrected t-test results
        fw_display_average_corr_threshmap = 0; % thresholded map corrected t-test results (Bonferroni)
        
        % individual maps and stats
        fw_display_all_zmaps = 0; % z-standardised average FWs
        fw_display_all_uncorr_thresh_maps = 0; % thresholded map uncorrected t-test results
        fw_display_all_corr_thresh_maps = 0; % thresholded map corrected t-test results (Bonferroni)
%__________________________________________________________________________    

elseif input_mode == 1 % Prompted manual input
    
    % specify analysis channels
    allchan = input('Are all possible channels analysed? "0" for no; "1" for yes (default if spatial/spatio-temporal): ');
    
    if allchan ~= 1
        
        relchan = input('Enter the channels to be analysed (e.g. [1 4 5]): ');
        
    end
    
    % specify properties of the decoding analysis
    stmode = input('Specify the s/t-analysis mode of the original analysis. "1" spatial, "2" temporal. "3" spatio-temporal: ');
    avmode = input('Specify the average mode of the original analysis. "1" single-trial, "2" run-average: ');
    window_width_ms = input('Specify the window width [ms] of the original analysis: ');
    step_width_ms = input('Specify the step width [ms] of the original analysis: ');
    
    % specify stats
    permstats = input('Testing against: "1" chance level; "2" permutation distribution: ');
    
    if permstats == 2
        
        drawmode = input('Testing against: "1" average permutated distribution (default); "2" random values drawn form permuted distribution (stricter): ');
        permdisp = input('Do you wish to display chance-level test results in figure? "0" for no; "1" for yes: ');
        
    end
    
    pstats = input('Specify critical p-value for statistical testing (e.g. 0.05): ');
    multcompstats = input(['\nSpecify if you wish to control for multiple comparisons: \n"0" for no correction \n'...
        '"1" for Bonferroni \n"2" for Holm-Bonferroni \n"3" for Strong FWER Control Permutation Testing \n' ...
        '"4" for Cluster-Based Permutation Testing \n"5" for KTMS Generalised FWER Control \n' ...
        '"6" for Benjamini-Hochberg FDR Control \n"7" for Benjamini-Krieger-Yekutieli FDR Control \n' ...
        '"8" for Benjamini-Yekutieli FDR Control \n Option: ']);
    
    if multcompstats == 3 || multcompstats == 4 || multcompstats == 5 % For permutation tests
        nIterations = input('Number of permutation iterations for multiple comparisons procedure (at least 1000 is recommended): ');    
    end
    if multcompstats == 5 % For KTMS Generalised FWER control
       KTMS_u = input('Enter the u parameter for the KTMS Generalised FWER control procedure: '); 
    end
    if multcompstats == 4 % For cluster-based permutation testing
       cluster_test_alpha = input('Enter the clustering threshold for detecting effects at individual time points (e.g. 0.05): '); 
    end
    
    % specify display options
    disp_on = input('Do you wish to display the results in figure(s)? "0" for no; "1" for yes: ');
    disp_sign = input('Specify if you wish to highlight significant results in figure. "0" for no; "1" for yes: ');
    
    % analyse feature weights
    fw_do = input('Do you wish to analyse the feature weights (only for spatial or spatio-temporal decoding)? "0" for no; "1" for yes: ');
    
    if fw_do == 1
        multcompstats = input(['\nSpecify which multiple comparisons correction method to use: \n' ...
        '"1" for Bonferroni \n"2" for Holm-Bonferroni \n"3" for Strong FWER Control Permutation Testing \n' ...
        '"4" for Cluster-Based Permutation Testing (Currently not available) \n"5" for KTMS Generalised FWER Control \n' ...
        '"6" for Benjamini-Hochberg FDR Control \n"7" for Benjamini-Krieger-Yekutieli FDR Control \n' ...
        '"8" for Benjamini-Yekutieli FDR Control \n Option: ']);
    
        if multcompstats == 3 || multcompstats == 4 || multcompstats == 5 % For permutation tests
            nIterations = input('Number of permutation iterations for multiple comparisons procedure (at least 1000 is recommended): ');    
        end
        if multcompstats == 5 % For KTMS Generalised FWER control
           KTMS_u = input('Enter the u parameter for the KTMS Generalised FWER control procedure: '); 
        end
        if multcompstats == 4 % For cluster-based permutation testing
           fprintf('Cluster-based corrections are currently not available.\n')
           % cluster_test_alpha = input('Enter the clustering threshold for detecting effects at individual time points (e.g. 0.05): '); 
        end
        
        fw_display_average_zmap = input('Do you wish to display the group-level averaged, z-standardised feature weights as a heat map? "0" for no; "1" for yes: '); % z-standardised average FWs
        fw_display_average_uncorr_threshmap = input(...
            'Do you wish to display the statistical threshold map (uncorrected) for the group-level averaged, z-standardised feature weights as a heat map? "0" for no; "1" for yes: '); % thresholded map uncorrected t-test results
        fw_display_average_corr_threshmap = input(...
            'Do you wish to display the statistical threshold map (corrected for multiple comparisons) for the group-level averaged, z-standardised feature weights as a heat map? "0" for no; "1" for yes: '); % thresholded map corrected t-test results (Bonferroni)
        
        % individual maps and stats
        fw_display_all_zmaps = input('');
        fw_display_all_uncorr_thresh_maps = input(...
            'Do you wish to display the statistical threshold map (uncorrected) for the group-level z-standardised feature weights for each time-step as a heat map? "0" for no; "1" for yes: ');
        fw_display_all_corr_thresh_maps = input(...
            'Do you wish to display the statistical threshold map (corrected for multiple comparisons) for the group-level z-standardised feature weights for each time-step as a heat map? "0" for no; "1" for yes: ');
        
    end
    
end % input

% Other inputs:
sbjs_todo = [1:5]; % Subject datasets to process
dcg_todo = 1; % Discrimination group results to analyse

% Call to function with all variadic inputs:
ANALYSE_DECODING_ERP(SLIST, sbjs_todo, dcg_todo, ...
    'allchan', allchan,...
    'relchan', relchan,...
    'stmode', stmode,...
    'avmode', avmode,...
    'window_width_ms', window_width_ms,...
    'step_width_ms', step_width_ms,...
    'permstats', permstats,...
    'drawmode', drawmode,...
    'pstats', pstats,...
    'multcompstats', multcompstats,...
    'nIterations', nIterations,...
    'KTMS_u', KTMS_u,...
    'cluster_test_alpha', cluster_test_alpha,...
    'disp_on', disp_on,...
    'permdisp', permdisp,...
    'disp_sign', disp_sign,...
    'fw_do', fw_do,...
    'fw_multcompstats', fw_multcompstats,...
    'fw_display_matrix', fw_display_matrix,...
    'fw_display_average_zmap', fw_display_average_zmap,...
    'fw_display_average_uncorr_threshmap', fw_display_average_uncorr_threshmap,...
    'fw_display_average_corr_threshmap', fw_display_average_corr_threshmap,...
    'fw_display_all_zmaps', fw_display_all_zmaps,...
    'fw_display_all_uncorr_thresh_maps', fw_display_all_uncorr_thresh_maps,...
    'fw_display_all_corr_thresh_maps', fw_display_all_corr_thresh_maps...
    );






