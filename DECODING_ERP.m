function DECODING_ERP(study_name,vconf,input_mode,sbj,dcg_todo,cross)
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
% This master-script prepares the data for the decoding analysis. It requires a
% configuration-script (e.g. "DEMO_config_v1.m") for the current study in 
% which all study parameters are specified. 
%
% The script outputs the average decoding accuracies for each analysis & time-
% step as well as feature weights.
%
% This toolbox interacts with LIBSVM toolbox (Chang & Lin) to do the classfication / regression
% see: https://www.csie.ntu.edu.tw/~cjlin/libsvm/
% Chang CC, Lin CJ (2011). LIBSVM : a library for support vector machines. ACM TIST, 2(3):27,
%
% requires:
% - study_name (e.g. 'DEMO')
% - vconf (version of study configuration script, e.g., "1" for DEMO_config_v1)
% - input_mode (1 = use coded varialbles from first section / 2 = enter manually)
% - sbj (number of subject to analyse, e.g., 1)
% - dcg_todo (discrimination group to analyse, as specified in SLIST.dcg_labels{dcg})
% - cross (cross-condition classification (0 = no / 1 = yes): Allows training on
% data from one DCG and predicting the left-out data from the respective
% condition in the other DCG (if "1" is chosen, two DCDs have to be entered
% for dcg_todo, e.g. [1 2])

%__________________________________________________________________________
%
% Variable naming convention: STRUCTURE_NAME.example_variable

%% SECTION 1: PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

global SBJTODO;
SBJTODO = sbj;
STUDY.sbj = sbj;

global DCGTODO;
DCGTODO = dcg_todo;
STUDY.dcg_todo = dcg_todo;
STUDY.cross = cross; % cross-condition decoding: 0 = no; 1 = A=>B; 2 = B=>A;

% decoding script = 2 - needed for interaction with other scripts to regulate
% functions such as saving data, calling specific sub-sets of parameters
global CALL_MODE
CALL_MODE = 2;

% get subject parameters
% sbj_list = input('Enter subject list: ','s');
STUDY.study_name = study_name;
STUDY.vconf = vconf;
STUDY.sbj_list = [study_name '_config_v' num2str(vconf)]; % use latest slist-function!

global SLIST;
eval(STUDY.sbj_list); % Gets subject parameters from the configuration script
%__________________________________________________________________________

% get study parameters and save in STUDY structure
%__________________________________________________________________________

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

%__________________________________________________________________________
% Adjust window and step widths using the sampling rate
STUDY.sampling_rate = SLIST.sampling_rate;
STUDY.window_width = floor(STUDY.window_width_ms / ((1/STUDY.sampling_rate) * 1000));
STUDY.step_width = floor(STUDY.step_width_ms / ((1/STUDY.sampling_rate) * 1000));

% cross-validation defaults for single-trial decoding
STUDY.n_all_analyses = STUDY.cross_val_steps * STUDY.n_rep_cross_val;
STUDY.n_all_permutation = STUDY.cross_val_steps * STUDY.permut_rep;

if STUDY.cross == 0 
    fprintf('DCG %d will be analysed. \n',STUDY.dcg_todo);                
elseif STUDY.cross > 0
    fprintf('DCG %d and %d will be analysed for cross-condition classification. \n',STUDY.dcg_todo(1), STUDY.dcg_todo(2));
end
    
fprintf('Cross-validation defaults for single-trial decoding:\n');
fprintf('%d steps with %d cycles resulting in %d analyses.\n',STUDY.cross_val_steps,STUDY.n_rep_cross_val,STUDY.n_all_analyses);
fprintf('Balanced number of examples will be used for all analyses by default.\n');

if STUDY.perm_test == 1
    fprintf('Random-label analysis will be based on %d analyses.\n',STUDY.n_all_permutation);
end


%% SECTION 2: READ IN DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic data in: eeg_sorted_cond{run,cond}(timepoints,channels,trials)
% *** converted into work_data{run,cond}(timepoints,channels,trials)

fprintf('Reading in data. Please wait... \n');
open_name = (SLIST.data_open_name);
load(open_name);
fprintf('Data loading complete.\n');

% read in regression labels if performing continuous regression
if STUDY.analysis_mode == 4 %i.e. if we are performing continuous regression
    
   STUDY.training_label_open_name = SLIST.training_label_file;
   STUDY.training_label_import = load(STUDY.training_label_open_name);
   STUDY.test_label_open_name = SLIST.test_label_file;
   STUDY.test_label_import = load(STUDY.test_label_open_name);
   fprintf('Loaded regressand labels.\n');

end

%__________________________________________________________________________
% if rt_match==1
%     
%     % match_my_rt.m is a toolbox-script that finds RT-matched trials in two
%     % conditions form dcg_todo. Only these conditions are matched, the
%     % other condition's data is carried forward but not used for decoding.
%     fprintf('Peforming RT-match for conditions.\n');
%     [eeg_sorted_cond_matched_rt,ntrials_dcg_todo,final_rts]=match_my_rts(eeg_sorted_cond,rt_sorted_cond,STUDY.dcg_todo);
%     clear eeg_sorted_cond;
%     eeg_sorted_cond=eeg_sorted_cond_matched_rt;
%     clear eeg_sorted_cond_matched_rt;
%   
% end
%__________________________________________________________________________

work_data = eval(SLIST.data_struct_name); % Copies eeg_sorted_cond data into work_data

%__________________________________________________________________________
% if data is stored in a structure and not a cell:
if isstruct(work_data)
    fprintf('Converting EEG-structure into cell.\n');
    wd = struct2cell(work_data); clear work_data;
    for i = 1:size(wd,2)
        for j = 1:size(wd,3)
            temp(:,:,:) = wd{1,i,j}(:,:,:);
            temp = double(temp);
            work_data{i,j} = temp;
            clear temp;
        end
    end
end
%__________________________________________________________________________

clear wd;
clear eeg_sorted_cond;

%% SECTION 3: REDUCE TO SPECIFIED DCG / conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STUDY_slist defines all DCG with their respective relevant conditions. In
% this section, work_data is reduced to the specified conditions
% *** work_data will be reduced to reduced_data{dcg,run,cond}(timepoints,channels,trials)
% (this section was modified in v29 to accomodate cross-condition decoding. A dimension was added to the data-matrix
%__________________________________________________________________________

if STUDY.analysis_mode ~= 4 % continuous SVR does not require conditions
          
    for r = 1:size(work_data,1) % for all runs
        
        % go through either one DCG (regular) or two DCGs (cross-condition decoding)
        for d = 1:size(STUDY.dcg_todo,2)
            
            for cond = 1:size(SLIST.dcg{STUDY.dcg_todo(d)},2) % for all conditions specified

                fprintf('Run %d: Extracting condition %d as specified in DCG %d.\n',r,(SLIST.dcg{STUDY.dcg_todo(d)}(cond)),STUDY.dcg_todo(d));
                temp(:,:,:) = work_data{r,(SLIST.dcg{STUDY.dcg_todo(d)}(cond))}(:,:,:);
                reduced_data{d,r,cond} = temp; 

                clear temp;                    

            end % cond    
        end % d
    end % r 

elseif STUDY.analysis_mode == 4 % put in by Dan 14/3/2014
    reduced_data{1,:,:} = work_data{:,:}; 
end % of if STUDY.analysis_mode ~= 4 statement

STUDY.nconds = size(reduced_data,3); % Calculates the number of conditions
clear work_data;

%% SECTION 4: DELETE EXTRA CHANNELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% these can be eye channels or remaining extra channels that have not been
% removed during pre-processing 
% *** shuffle_data will be converted to clean_data{dcg,run,cond}(timepoints,channels,trials)
% (In v29 a dimension was added to accomodate cross-condition decoding)
%__________________________________________________________________________

for d = 1:size(reduced_data,1)
    
    for r = 1:size(reduced_data,2)

        for cond = 1:size(reduced_data,3)

            temp_data = reduced_data{d,r,cond};
            szt = size(temp_data); % DF NOTE: This variable szt appears to be unused.
            STUDY.nextra = SLIST.extra;

            if STUDY.nextra > 0
                fprintf('Removing extra channels for participant %d condition %d \n',SBJTODO,cond);
                extra_channels = SLIST.extra;
                temp_data(:,extra_channels,:) = [];
            elseif STUDY.nextra == 0      
                fprintf('No extra channels removed for participant %d, run %d, condition %d \n',SBJTODO,r,cond);
            end

            clean_data{d,r,cond} = temp_data;
            clear temp_data;     
        end % cond
    end % run
end %d
clear reduced_data;

%% SECTION 5: CALCULATE MIN NUMBER TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find minimum number of trials for conditions in discrimination group
% will be the same in RT-matched data anyway, otherwise random selection due to shuffling
% *** clean_data will be converted to balanced_data{dcg,run,cond}(timepoints,channels,trials)
% (In v29 a dimension was added to accomodate cross-condition decoding)
%__________________________________________________________________________

% calculate min number of trials
for r = 1:size(clean_data,2)
    
    mintrs = [];
    for d = 1:size(clean_data,1)
        
        for cond = 1:size(clean_data,3)

            temp = clean_data{d,r,cond};
            ntrials_cond = size(temp,3);
            mintrs = [mintrs ntrials_cond];
            clear temp;

        end % cond
       
    end % d 
    
    % minimum per run (all DCGs)
    STUDY.mintrs_min(r) = min(mintrs);
    
end %r

fprintf('Minimum number of trials per condition computed for participant %d \n',SBJTODO);

% select min trials in all conditions
for r = 1:size(clean_data,2)
    
    for d = 1:size(clean_data,1)

        for cond = 1:size(clean_data,3)

            % extract min number trials for this DCG and run
            temp_data = clean_data{d,r,cond};
            balanced_data{d,r,cond} = temp_data(:,:,1:(STUDY.mintrs_min(r)));
            clear temp_data;

        end % cond          

    end % r
end %d
fprintf('Same number of trials is used for all conditions within each run.\n');

clear clean_data;

%% SECTION 6: CALCULATE RUN-AVERAGES / POOL DATA ACROSS RUNS %%%%%%%%%%%%%%%%%%%%%%%%%
% data is either averaged across trials within runs: mean_balanced_data{run,cond}(timepoints,channels){run,condition}
% or pooled across runs (if available): pooled_balanced_data{DCG,1,cond}(timepoints,channels,trials)
% if only one run available, data will be unchnaged: balanced_data will be
% replaced by mean_balanced_data
% (In v29 a dimension was added to accomodate cross-condition decoding)
%__________________________________________________________________________

if STUDY.avmode == 1 % single-trials 
    
    for d = 1:size(balanced_data,1)
        for cond = 1:size(balanced_data,3)
            for r = 1:size(balanced_data,2)

                if r == 1
                    all_data = balanced_data{d,r,cond};
                elseif r > 1
                    all_data = cat(3,all_data,balanced_data{d,r,cond});
                end

            end % r
            pooled_balanced_data{d,1,cond} = all_data;
            clear all_data;
        end % cond
    end %d
    fprintf('Data from all runs (if more than one) have been pooled into one dataset.\n');    
    
elseif STUDY.avmode == 2 % run averages
    
    for d = 1:size(balanced_data,1)
        for r = 1:size(balanced_data,2)
            for cond = 1:size(balanced_data,3)

                mean_balanced_data{d,r,cond} = mean(balanced_data{d,r,cond},3);

            end % cond        
        end % r
    end %d
    fprintf('Run averages based across trials were computed for each condition.\n');
end % of if STUDY.avmode == 1 statement

clear balanced_data;

%% SECTION 7: SORT DATA FOR CLASSIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output of this section:
%
% train_set{dcg, condition, cross-validation step, cycles of cross-validation}(datapoints, channels, trials/exemplar)
% test_set{dcg, condition, cross-validation step, cycles of cross-validation}(datapoints, channels, trials/exemplar)
%
% for run-averaged data the cross-validation steps = number of runs and
% one cicle of cross-validation only
% every cell contains one data-points x channels (x exemplars for one classification 
% step = trials / not necessary for run-average decoding) matrix
% (In v29 a dimension was added to accomodate cross-condition decoding)
%__________________________________________________________________________

% single-trial data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all runs are pooled. Data is randomly drawn from all data
if STUDY.avmode == 1
   
    for d = 1:size(pooled_balanced_data,1)
        
        % take X% of data as test data with X-fold cross-validation
        % repeat X-times with different test data sets 
        for con = 1:size(pooled_balanced_data,3)

            ntrs_set = floor(size(pooled_balanced_data{d,1,con}(:,:,:),3) / STUDY.cross_val_steps);

            for ncv = 1:STUDY.n_rep_cross_val

                % draw x times without replacement
                random_set = randperm(ntrs_set * STUDY.cross_val_steps);
                x = 1;

                for cv = 1:STUDY.cross_val_steps

                    % find the test-trials for current cross-validation step
                    test_trials = random_set(x:(x + (ntrs_set - 1)));

                    temp_training_set = pooled_balanced_data{d,1,con}(:,:,:);

                    if STUDY.analysis_mode == 4 % if continuous SVR
                        temp_training_labels = STUDY.training_label_import.labels{1,con};
                    end % if

                    delete_test_trials = fliplr(sort(test_trials));

                    % extract test-trials and delete them from training-trials
                    for trl = 1:size(test_trials,2)

                        test_set{d,con,cv,ncv}(:,:,trl) = pooled_balanced_data{d,1,con}(:,:,test_trials(trl));

                        if STUDY.analysis_mode == 4 % if continuous SVR
                            STUDY.test_labels{con,cv,ncv}(trl) = STUDY.test_label_import.labels{1,con}(test_trials(trl));
                        end % if

                        % option 1: delete used trials each for trial
                        % temp_training_set(:,:,delete_test_trials(trl))=[];

                    end %trl

                    % option 2: delete all used trials at once
                    temp_training_set(:,:,delete_test_trials) = [];

                    training_set{d,con,cv,ncv} = temp_training_set;

                    if STUDY.analysis_mode == 4 % if continuous SVR
                        temp_training_labels(delete_test_trials) = [];
                        STUDY.training_labels{con,cv,ncv} = temp_training_labels;
                    end

                    clear delete_test_trials;
                    clear temp_training_set;
                    clear temp_training_labels;
                    clear test_trials;

                    % go to next step and repeat
                    x = x + ntrs_set;
                    fprintf('Data sorted for single-trial decoding: specified DCG %d condition %d cross-validation step %d of cycle %d. \n',d,con,cv,ncv);    

                end % cv (cross-validation steps)

                clear random_set

            end % ncv (repetition of cross-validation)

            clear ntrs_set;

        end % con (all conditions in serial order)
    
    end % d
    clear pooled_balanced_data
    
% run-average data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% test- and training-data correspond to run averages 
elseif STUDY.avmode == 2 % run averages
    
    for d = 1:size(mean_balanced_data,1)
        
        for con = 1:size(mean_balanced_data,3)

            % take data from each run as test data once with number of runs =
            % cross-validation steps. Don't repeat because data is not randomly drawn from experiment        
            for r = 1:size(mean_balanced_data,2)

                % r = test data-set
                test_run = r;   % DF NOTE: Variable test_run appears to be unused       
                test_set{d,con,r,1} = mean_balanced_data{d,r,con};

                % all other runs = training data-set
                train_on = find(1:(size(mean_balanced_data,2)) ~= r);
                sz_train = size(train_on,2);

                for trainrun = 1:sz_train

                    if trainrun == 1
                        train_data = mean_balanced_data{d,train_on(trainrun),con};
                    elseif trainrun > 1     
                        train_data = cat(3,train_data,mean_balanced_data{d,train_on(trainrun),con});   
                    end

                end % trainrun

                training_set{d,con,r,1} = train_data;
                clear train_data;
                fprintf('Data sorted for run-average decoding: condition %d cross-validation step %d. \n',con,r);   

            end % r

        end % con
    
    end % d
    
    clear mean_balanced_data
    
end % avmode

%% SECTION 8: BUILD LABELS AND DO CLASSIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

% input('Do you wish to continue with decoding? Press any button!')

fprintf('Starting with vector preparation... \n');
[RESULTS] = prepare_my_vectors_erp(training_set, test_set, SLIST, STUDY);

%% SECTION 9: AVERAGE CROSS-VALIDATION STEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

RESULTS.subj_acc = [];
RESULTS.subj_perm_acc = [];

% na = number analyses (stmode=2 has one per channel)
for na = 1:size(RESULTS.prediction_accuracy,2)
    
    pa(:,:,:) = RESULTS.prediction_accuracy{na}(:,:,:);
    
    if STUDY.perm_test == 1
        perm_pa(:,:,:) = RESULTS.perm_prediction_accuracy{na}(:,:,:);
    end
               
    % calculate average decoding accuracy
    RESULTS.subj_acc(na,:) = nanmean(nanmean(pa,3),2);
    clear pa;
        
    % calculate average permutation test decoding accuracy
    if STUDY.perm_test == 1
        RESULTS.subj_perm_acc(na,:) = nanmean(nanmean(perm_pa,3),2);
        clear perm_pa;
    end
        
end % na

fprintf('Results are computed and averaged for participant %d. \n',STUDY.sbj);

%% SECTION 10: SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saves decoding results to a .mat file in the output directory. Some analysis
% settings (e.g. window_width_ms) are included in the file name.
%__________________________________________________________________________

if STUDY.cross == 0
    savename = [(SLIST.output_dir) STUDY.study_name '_SBJ' num2str(STUDY.sbj) '_win' num2str(STUDY.window_width_ms) '_steps' num2str(STUDY.step_width_ms)...
        '_av' num2str(STUDY.avmode) '_st' num2str(STUDY.stmode) '_DCG' SLIST.dcg_labels{STUDY.dcg_todo} '.mat'];
elseif STUDY.cross == 1
    savename = [(SLIST.output_dir) STUDY.study_name '_SBJ' num2str(STUDY.sbj) '_win' num2str(STUDY.window_width_ms) '_steps' num2str(STUDY.step_width_ms)...
        '_av' num2str(STUDY.avmode) '_st' num2str(STUDY.stmode) '_DCG' SLIST.dcg_labels{STUDY.dcg_todo(1)} 'toDCG' SLIST.dcg_labels{STUDY.dcg_todo(2)} '.mat'];
elseif STUDY.cross == 2
        savename = [(SLIST.output_dir) STUDY.study_name '_SBJ' num2str(STUDY.sbj) '_win' num2str(STUDY.window_width_ms) '_steps' num2str(STUDY.step_width_ms)...
        '_av' num2str(STUDY.avmode) '_st' num2str(STUDY.stmode) '_DCG' SLIST.dcg_labels{STUDY.dcg_todo(2)} 'toDCG' SLIST.dcg_labels{STUDY.dcg_todo(1)} '.mat'];
end
    
save(savename,'STUDY','RESULTS'); % Save STUDY and RESULTS structures into a .mat file

fprintf('Results are saved for participant %d in directory: %s. \n',STUDY.sbj,(SLIST.output_dir));

%% SECTION 11: DISPLAY INDIVIDUAL RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displays decoding results for each individual subject dataset.
%__________________________________________________________________________

if STUDY.display_on == 1
    
    display_indiv_results_erp(STUDY,RESULTS);

end % display_on
%__________________________________________________________________________ 