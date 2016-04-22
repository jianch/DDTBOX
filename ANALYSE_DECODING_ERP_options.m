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
    ANALYSIS.multcompstats = 0; % Correction for multiple comparisons: 
    % 0 = no correction
    % 1 = Bonferroni correction
    % 2 = Holm-Bonferroni correction
    % 3 = Strong FWER Control Permutation Test
    % 4 = Cluster-Based Permutation Test
    % 5 = KTMS Generalised FWER Control
    % 6 = Benjamini-Hochberg FDR Control
    % 7 = Benjamini-Krieger-Yekutieli FDR Control
    % 8 = Benjamini-Yekutieli FDR Control
    ANALYSIS.nIterations = 5000; % Number of permutation or bootstrap iterations for resampling-based multiple comparisons correction procedures
    ANALYSIS.KTMS_u = 2; % u parameter of the KTMS GFWER control procedure
    ANALYSIS.cluster_test_alpha = 0.05; % For cluster-based test: Significance threshold for detecting effects at individual time windows (e.g. 0.05)
    
    ANALYSIS.disp.on = 1; % display a results figure? 0=no / 1=yes
    ANALYSIS.permdisp = 1; % display the results from permutation test in figure as separate line? 0=no / 1=yes
    ANALYSIS.disp.sign = 1; % display statistically significant steps in results figure? 0=no / 1=yes
    
    ANALYSIS.fw.do = 0; % analyse feature weights? 0=no / 1=yes
    ANALYSIS.fw.multcompstats = 1; % Feature weights correction for multiple comparisons:
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
        ANALYSIS.fw.display_matrix = 1; % feature weights matrix % 0=no / 1=yes
        
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
    
    if ANALYSIS.multcompstats == 3 || ANALYSIS.multcompstats == 4 || ANALYSIS.multcompstats == 5 % For permutation tests
        ANALYSIS.nIterations = input('Number of permutation iterations for multiple comparisons procedure (at least 1000 is recommended): ');    
    end
    if ANALYSIS.multcompstats == 5 % For KTMS Generalised FWER control
       ANALYSIS.KTMS_u = input('Enter the u parameter for the KTMS Generalised FWER control procedure: '); 
    end
    if ANALYSIS.multcompstats == 4 % For cluster-based permutation testing
       ANALYSIS.cluster_test_alpha = input('Enter the clustering threshold for detecting effects at individual time points (e.g. 0.05): '); 
    end
    
    % specify display options
    ANALYSIS.disp.on = input('Do you wish to display the results in figure(s)? "0" for no; "1" for yes: ');
    ANALYSIS.disp.sign = input('Specify if you wish to highlight significant results in figure. "0" for no; "1" for yes: ');
    
    % analyse feature weights
    ANALYSIS.fw.do = input('Do you wish to analyse the feature weights (only for spatial or spatio-temporal decoding)? "0" for no; "1" for yes: ');
    
    if ANALYSIS.fw.do == 1
        ANALYSIS.multcompstats = input(['\nSpecify which multiple comparisons correction method to use: \n' ...
        '"1" for Bonferroni \n"2" for Holm-Bonferroni \n"3" for Strong FWER Control Permutation Testing \n' ...
        '"4" for Cluster-Based Permutation Testing (Currently not available) \n"5" for KTMS Generalised FWER Control \n' ...
        '"6" for Benjamini-Hochberg FDR Control \n"7" for Benjamini-Krieger-Yekutieli FDR Control \n' ...
        '"8" for Benjamini-Yekutieli FDR Control \n Option: ']);
    
        if ANALYSIS.multcompstats == 3 || ANALYSIS.multcompstats == 4 || ANALYSIS.multcompstats == 5 % For permutation tests
            ANALYSIS.nIterations = input('Number of permutation iterations for multiple comparisons procedure (at least 1000 is recommended): ');    
        end
        if ANALYSIS.multcompstats == 5 % For KTMS Generalised FWER control
           ANALYSIS.KTMS_u = input('Enter the u parameter for the KTMS Generalised FWER control procedure: '); 
        end
        if ANALYSIS.multcompstats == 4 % For cluster-based permutation testing
           fprintf('Cluster-based corrections are currently not available.\n')
           % ANALYSIS.cluster_test_alpha = input('Enter the clustering threshold for detecting effects at individual time points (e.g. 0.05): '); 
        end
        
        ANALYSIS.fw.display_average_zmap = input('Do you wish to display the group-level averaged, z-standardised feature weights as a heat map? "0" for no; "1" for yes: '); % z-standardised average FWs
        ANALYSIS.fw.display_average_uncorr_threshmap = input(...
            'Do you wish to display the statistical threshold map (uncorrected) for the group-level averaged, z-standardised feature weights as a heat map? "0" for no; "1" for yes: '); % thresholded map uncorrected t-test results
        ANALYSIS.fw.display_average_corr_threshmap = input(...
            'Do you wish to display the statistical threshold map (corrected for multiple comparisons) for the group-level averaged, z-standardised feature weights as a heat map? "0" for no; "1" for yes: '); % thresholded map corrected t-test results (Bonferroni)
        
        % individual maps and stats
        ANALYSIS.fw.display_all_zmaps = input('');
        ANALYSIS.fw.display_all_uncorr_thresh_maps = input(...
            'Do you wish to display the statistical threshold map (uncorrected) for the group-level z-standardised feature weights for each time-step as a heat map? "0" for no; "1" for yes: ');
        ANALYSIS.fw.display_all_corr_thresh_maps = input(...
            'Do you wish to display the statistical threshold map (corrected for multiple comparisons) for the group-level z-standardised feature weights for each time-step as a heat map? "0" for no; "1" for yes: ');
        
    end
    
end % input