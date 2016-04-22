if input_mode == 0 % Hard-coded input

    % get study parameters (hardcoded)
    STUDY.stmode = 3; % SPACETIME mode (1=spatial / 2=temporal / 3=spatio-temporal)
    STUDY.avmode = 1; % AVERAGE mode (1=no averaging; single-trial / 2=run average) 
    STUDY.analysis_mode = 1;% ANALYSIS mode (1=SVM classification, 2=LDA classification, 3=SVR increments, 4=SVR continous)
    
    STUDY.window_width_ms = 10; % width of sliding window in ms
    STUDY.step_width_ms = 10; % step size with which sliding window is moved through the trial
    
    STUDY.perm_test = 1; % run the permutation-decoding? 0=no / 1=yes
    STUDY.perm_disp = 1; % display the permutation results in figure? 0=no / 1=yes
    STUDY.display_on = 1; % 1=figure displayed, 0=no figure
    
    STUDY.rt_match = 0; % Use RT-matching algorithm to select trials? 0=no / 1=yes
    STUDY.feat_weights_mode = 1; % Extract feature weights? 0=no / 1=yes
    STUDY.cross_val_steps = 10; % How many cross-validation steps (if no runs available)?
    STUDY.n_rep_cross_val = 10; % How many repetitions of full cross-validation with randomly re-ordered data?
    STUDY.permut_rep = 10; % How many repetitions of full cross-validation with permutation results?

elseif input_mode == 1 % Prompt user for input

    % get study parameters via input
    STUDY.stmode = input('Enter "1" for spatial, "2" for temporal, or "3" for spatio-temporal decoding: ');
    STUDY.avmode = input('Enter "1" for single-trial decoding or "2" for run-average decoding: ');
    if STUDY.avmode == 1 % Options for single-trial decoding
        STUDY.cross_val_steps = input('How many cross-validation steps do you wish to perform?');
        STUDY.n_rep_cross_val = input('How many independent repetitions of the analysis do you wish to perform?');
        STUDY.rt_match = input('Do you wish to RT-match the trials from the conditions to be decoded? (1=yes, 0=no) ');
    end
    STUDY.analysis_mode = input('Specifiy analysis method: "1" for Class SVM, "2" for Class LDA, "3" increments SVR, "4" continuous SVR: '); 
    STUDY.window_width_ms = input('Enter decoding window width in ms: ');
    STUDY.step_width_ms = input('Enter step width for moving the decoding window in ms: ');
    STUDY.cross = input('Do you wish to perform cross-condition decoding? "0" for no, "1" for yes:');
    if STUDY.cross > 0
        dcgs = input('Enter two discriminations groups for cross-decoding (e.g.[1 2]):');
        STUDY.dcg_todo = dcgs;
    end
        STUDY.perm_test = input('Do you want to run a permutation test? (1=yes, 0=no) '); 
    if STUDY.perm_test == 1
        STUDY.permut_rep = input('How many repetitions of the permutation test do you wish to perform? ');
        STUDY.perm_disp = input('Do you wish to plot the permutation test results? (1=yes, 0=no) '); 
    end
    STUDY.display_on = input('Do you wish to plot the individual decoding results? (1=yes, 0=no) ');

end