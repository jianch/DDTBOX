function decoding_erp(cfg)
%
% Performs MVPA on a single subject dataset.
%
%
% Inputs:
%
%   cfg structure containing subject dataset information and decoding analysis
%   parameters. Information about each parameter is described in the
%   project wiki and in the example configuration script 
%   (e.g. EXAMPLE_run_decoding_analyses.m)
%
%   
%		
% Optional keyword inputs:
%
%
%
% Usage:           decoding_erp(cfg)
%
%
% Copyright (c) 2013-2016, Stefan Bode and contributors 
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



%% SECTION 1: PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

cfg.sbj_todo = cfg.sbj;
STUDY.sbj = cfg.sbj;
STUDY.dcg_todo = cfg.dcg_todo;
STUDY.regr_todo = [];
STUDY.cross = cfg.cross; % cross-condition decoding: 0 = no; 1 = A=>B; 2 = B=>A;

% get subject parameters
% sbj_list = input('Enter subject list: ','s');
STUDY.study_name = cfg.study_name;

%__________________________________________________________________________

% get study parameters and save in STUDY structure
%__________________________________________________________________________

% List of parameters for updating header 
% TODO: Update header with all cfg
% parameters.
%     % Multivariate classification/regression parameters
%     STUDY.analysis_mode = cfg.analysis_mode; % ANALYSIS mode (1=SVM with LIBSVM / 2=SVM with liblinear / 3=SVR with LIBSVM)
%     STUDY.stmode = cfg.stmode; % SPACETIME mode (1=spatial / 2=temporal / 3=spatio-temporal)
%     STUDY.avmode = cfg.avmode; % AVERAGE mode (1=no averaging; single-trial / 2=run average) !!! need single trials for SVR !!!
%     STUDY.window_width_ms = cfg.window_width_ms; % width of sliding window in ms
%     STUDY.step_width_ms = cfg.step_width_ms; % step size with which sliding window is moved through the trial
%     STUDY.zscore_convert = cfg.zscore_convert; % Convert data into z-scores before decoding? 0 = no / 1 = yes
%     STUDY.perm_test = cfg.perm_test; % run the permutation-decoding? 0=no / 1=yes
%     STUDY.cross_val_steps =  cfg.cross_val_steps; % How many cross-validation steps (if no runs available)?
%     STUDY.n_rep_cross_val = cfg.n_rep_cross_val; % How many repetitions of full cross-validation with re-ordered data?
%     STUDY.permut_rep = cfg.permut_rep; % How many repetitions of full cross-validation with permutation results?
%     
%     % Feature weights extraction
%     STUDY.feat_weights_mode = cfg.feat_weights_mode; % Extract feature weights? 0=no / 1=yes
% 
%     % Single subject decoding results plotting
%     STUDY.display_on = cfg.display_on; % Display individual subject results? 1=figure displayed / 0=no figure
%     STUDY.perm_disp = cfg.perm_disp; % display the permutation results in figure? 0=no / 1=yes

if cfg.analysis_mode == 1 
    cfg.analysis_mode_label = 'SVM_LIBSVM';
elseif cfg.analysis_mode == 2 
    cfg.analysis_mode_label = 'SVM_LIBLIN';
elseif cfg.analysis_mode == 3
    cfg.analysis_mode_label = 'SVR_LIBSVM';
end

%__________________________________________________________________________
% Adjust window and step widths using the sampling rate
cfg.sampling_rate = cfg.sampling_rate;
cfg.pointzero = cfg.pointzero;
cfg.dcg_label = cfg.dcg_labels{cfg.dcg_todo};

cfg.window_width = floor(cfg.window_width_ms / ((1/cfg.sampling_rate) * 1000));
cfg.step_width = floor(cfg.step_width_ms / ((1/cfg.sampling_rate) * 1000));

% cross-validation defaults for single-trial decoding
cfg.n_all_analyses = cfg.cross_val_steps * cfg.n_rep_cross_val;
cfg.n_all_permutation = cfg.cross_val_steps * cfg.permut_rep;

if cfg.cross == 0 
    fprintf('DCG %d will be analysed. \n',cfg.dcg_todo);                
elseif cfg.cross > 0
    fprintf('DCG %d and %d will be analysed for cross-condition classification. \n', cfg.dcg_todo(1), cfg.dcg_todo(2));
end

% For SVR dcg_todo defines condition that the target variable comes from 
% (to be saved in SLIST.regress_struct_name(rows=trials; columns = varialbe values)
% The data is always the same, so by default: dcg_todo = 1 
if cfg.analysis_mode == 3 
    cfg.regr_todo = cfg.dcg_todo;
    cfg.dcg_todo = 1;
end
    
fprintf('Cross-validation defaults for single-trial decoding:\n');
fprintf('%d steps with %d cycles resulting in %d analyses.\n',cfg.cross_val_steps,cfg.n_rep_cross_val,cfg.n_all_analyses);
fprintf('Balanced number of examples will be used for all analyses by default.\n');

if cfg.perm_test == 1
    fprintf('Random-label analysis will be based on %d analyses.\n',cfg.n_all_permutation);
end

%__________________________________________________________________________
% LIBSVM and LIBLINEAR flags
% LIBSVM and LIBLINEAR require input flags (strings) to specify the type of model that will be used for decoding.
% WARNING: Do not change these flags unless you really know what you are doing!!!

% For the full list of flags and their options see
% LIBSVM:    https://www.csie.ntu.edu.tw/~cjlin/libsvm/
% LIBLINEAR: https://www.csie.ntu.edu.tw/~cjlin/liblinear/
% The following flags are currently used as defaults in DDTBox

% -c cost : cost parameter
% The following backends are

% LIBSVM
% -s svm_type : set the type of support vector machine
% 	0 -- C-Support Vector Classification
% 	1 -- nu-Support Vector Classification
% 	2 -- one-class Support Vector Machine
% 	3 -- epsilon-Support Vector Regression
% 	4 -- nu-Support Vector Regression
% -t kernel_type : set type of kernel function
% 	0 -- linear: u'*v
% 	1 -- polynomial: (gamma*u'*v + coef0)^degree
% 	2 -- radial basis function: exp(-gamma*|u-v|^2)
% 	3 -- sigmoid: tanh(gamma*u'*v + coef0)

% LIBLINEAR
% -s svm_type:
%	 0 -- L2-regularized logistic regression (primal)
%	 1 -- L2-regularized L2-loss support vector classification (dual)
%	 2 -- L2-regularized L2-loss support vector classification (primal)
%	 3 -- L2-regularized L1-loss support vector classification (dual)
%	 4 -- support vector classification by Crammer and Singer
%	 5 -- L1-regularized L2-loss support vector classification
%	 6 -- L1-regularized logistic regression
%	 7 -- L2-regularized logistic regression (dual)
%	11 -- L2-regularized L2-loss support vector regression (primal)
%	12 -- L2-regularized L2-loss support vector regression (dual)
%	13 -- L2-regularized L1-loss support vector regression (dual)


% Defaults:
% Support Vector Classification with libsvm - '-s 0 -t 0 -c 1'
% Support Vector Regression with libsvm - '-s 3 -t 0 -c 0.1'
% Support Vector Regression (continuous) with libsvm - '-s 3 -t 0 -c 0.1'
% Support Vector Classification with liblinear - '-s 2 -c 1'

if cfg.analysis_mode == 1 % SVM classification with libsvm
    cfg.backend_flags.svm_type = 0;
    cfg.backend_flags.kernel_type = 0;
    cfg.backend_flags.cost = 1;
    cfg.backend_flags.extra_flags = []; % To input extra flag types
elseif cfg.analysis_mode == 2 % SVM with liblinear
    cfg.backend_flags.svm_type = 2;
    cfg.backend_flags.kernel_type = -1; % not valid for liblinear
    cfg.backend_flags.cost = 1;
    cfg.backend_flags.extra_flags = [];
elseif cfg.analysis_mode == 3 % SVR (regression)  with libsvm
    cfg.backend_flags.svm_type = 3;
    cfg.backend_flags.kernel_type = 0;
    cfg.backend_flags.cost = 0.1;
    cfg.backend_flags.extra_flags = []; 
end

% Merging all flags into a single string
cfg.backend_flags.all_flags = ['-s ' int2str(cfg.backend_flags.svm_type) ' -c ' num2str(cfg.backend_flags.cost)];

% LIBSVM specific options
if cfg.backend_flags.kernel_type ~= -1
	cfg.backend_flags.all_flags = [cfg.backend_flags.all_flags ' -t ' int2str(cfg.backend_flags.kernel_type)];
end

cfg.backend_flags.all_flags = [cfg.backend_flags.all_flags ' ' cfg.backend_flags.extra_flags];


%% SECTION 2: READ IN DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic data in: eeg_sorted_cond{run,cond}(timepoints,channels,trials)
% *** converted into work_data{run,cond}(timepoints,channels,trials)

fprintf('Reading in data. Please wait... \n');
open_name = (cfg.data_open_name);
load(open_name);
fprintf('Data loading complete.\n');

% read in regression labels if performing SVR
if cfg.analysis_mode == 3
    
    cfg.regress_open_name = cfg.regress_label_name;
    cfg.regress_data = load(cfg.regress_open_name);    

    fprintf('Loaded regressand labels.\n');

end

work_data = eval(cfg.data_struct_name); % Copies eeg_sorted_cond data into work_data
cfg.nchannels=cfg.nchannels;

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

%__________________________________________________________________________
% if data dimensions are flipped (e.g, when taken directly from EEGlab):

if size(work_data{1,1},1)==cfg.nchannels && size(work_data{1,1},2)~=cfg.nchannels  
    
    for r = 1:size(work_data,1)
        for c = 1:size(work_data,2)
            
            temp(:,:,:) = work_data{r,c}(:,:,:);
            temp = permute(temp,[2 1 3]);
            work_data{r,c} = temp;
            clear temp;  
            
        end % c        
    end % r
    fprintf('Data was converted into the correct format: eeg_sorted_cond{run,cond}(timepoints,channels,trials).\n');   

elseif size(work_data{1,1},1)~=cfg.nchannels && size(work_data{1,1},2)==cfg.nchannels  
    
    fprintf('Data seems to be in the correct format: eeg_sorted_cond{run,cond}(timepoints,channels,trials).\n');
    
elseif size(work_data{1,1},1)==cfg.nchannels && size(work_data{1,1},2)==cfg.nchannels
    
    fprintf('Number of channels = number of time points? Check whether data is in the correct format.\n');
    
elseif size(work_data{1,1},1)~=cfg.nchannels && size(work_data{1,1},2)~=cfg.nchannels 
    
    fprintf('Data might not be in the required format: eeg_sorted_cond{run,cond}(timepoints,channels,trials). \n');
    
end;
%__________________________________________________________________________


%% SECTION 3: REDUCE TO SPECIFIED DCG / conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STUDY_slist defines all DCG with their respective relevant conditions. In
% this section, work_data is reduced to the specified conditions
% *** work_data will be reduced to reduced_data{dcg,run,cond}(timepoints,channels,trials)
% (this section was modified in v29 to accomodate cross-condition decoding. A dimension was added to the data-matrix
%__________________________________________________________________________

if cfg.analysis_mode ~= 3 % SVR does not require conditions
          
    for r = 1:size(work_data,1) % for all runs
        
        % go through either one DCG (regular) or two DCGs (cross-condition decoding)
        for d = 1:size(cfg.dcg_todo,2)
            
            for cond = 1:size(cfg.dcg{cfg.dcg_todo(d)},2) % for all conditions specified

                fprintf('Run %d: Extracting condition %d as specified in DCG %d.\n',r,(cfg.dcg{cfg.dcg_todo(d)}(cond)),cfg.dcg_todo(d));
                temp(:,:,:) = work_data{r,(cfg.dcg{cfg.dcg_todo(d)}(cond))}(:,:,:);
                reduced_data{d,r,cond} = temp; 

                clear temp;                    

            end % cond    
        end % d
    end % r 

elseif cfg.analysis_mode == 3 % put in by Dan 14/3/2014
    reduced_data{1,:,:} = work_data{:,:}; 
end % of if STUDY.analysis_mode ~= 4 statement

cfg.nconds = size(reduced_data,3); % Calculates the number of conditions
clear work_data;


%% SECTION 4: CALCULATE MIN NUMBER TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find minimum number of trials for conditions in discrimination group
% will be the same in RT-matched data anyway, otherwise random selection due to shuffling
% *** clean_data will be converted to balanced_data{dcg,run,cond}(timepoints,channels,trials)
% (In v29 a dimension was added to accomodate cross-condition decoding)
%__________________________________________________________________________

% calculate min number of trials
for r = 1:size(reduced_data,2)
    
    mintrs = [];
    for d = 1:size(reduced_data,1)
        
        for cond = 1:size(reduced_data,3)

            temp = reduced_data{d,r,cond};
            ntrials_cond = size(temp,3);
            mintrs = [mintrs ntrials_cond];
            clear temp;

        end % cond
       
    end % d 
    
    % minimum per run (all DCGs)
    cfg.mintrs_min(r) = min(mintrs);
    
end %r

fprintf('Minimum number of trials per condition computed for participant %d \n', cfg.sbj_todo);

% select min trials in all conditions
for r = 1:size(reduced_data,2)
    
    for d = 1:size(reduced_data,1)

        for cond = 1:size(reduced_data,3)

            % extract min number trials for this DCG and run
            temp_data = reduced_data{d,r,cond};
            balanced_data{d,r,cond} = temp_data(:,:,1:(cfg.mintrs_min(r)));
            clear temp_data;

        end % cond          

    end % r
end %d
fprintf('Same number of trials is used for all conditions within each run.\n');

clear reduced_data;


%% SECTION 5: CALCULATE RUN-AVERAGES / POOL DATA ACROSS RUNS %%%%%%%%%%%%%%%%%%%%%%%%%
% data is either averaged across trials within runs: mean_balanced_data{run,cond}(timepoints,channels){run,condition}
% or pooled across runs (if available): pooled_balanced_data{DCG,1,cond}(timepoints,channels,trials)
% if only one run available, data will be unchnaged: balanced_data will be
% replaced by mean_balanced_data
% (In v29 a dimension was added to accomodate cross-condition decoding)
%__________________________________________________________________________

if cfg.avmode == 1 % single-trials 
    
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
    
elseif cfg.avmode == 2 % run averages
    
    for d = 1:size(balanced_data,1)
        for r = 1:size(balanced_data,2)
            for cond = 1:size(balanced_data,3)

                mean_balanced_data{d,r,cond} = mean(balanced_data{d,r,cond},3);

            end % cond        
        end % r
    end %d
    fprintf('Run averages based across trials were computed for each condition.\n');
end % of if cfg.avmode == 1 statement

clear balanced_data;


%% SECTION 6: SORT DATA FOR CLASSIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output of this section:
%
% training_set{dcg, condition, cross-validation step, cycles of cross-validation}(datapoints, channels, trials/exemplar)
% test_set{dcg, condition, cross-validation step, cycles of cross-validation}(datapoints, channels, trials/exemplar)
%
% for run-averaged data the cross-validation steps = number of runs and
% one cycle of cross-validation only
% every cell contains one data-points x channels (x exemplars for one classification 
% step = trials / not necessary for run-average decoding) matrix
% (In v29 a dimension was added to accomodate cross-condition decoding)
%__________________________________________________________________________

% single-trial data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all runs are pooled. Data is randomly drawn from all data
if cfg.avmode == 1
    
    % take X% of data as test data with X-fold cross-validation
    % repeat X-times with different test data sets
    for con = 1:size(pooled_balanced_data,3)
        
        ntrs_set = floor(size(pooled_balanced_data{d,1,con}(:,:,:),3) / cfg.cross_val_steps);
        
        for ncv = 1:cfg.n_rep_cross_val
            
            % draw x times without replacement
            random_set = randperm(ntrs_set * cfg.cross_val_steps);
            x = 1;
            
            for cv = 1:cfg.cross_val_steps
                
                % find the test-trials for current cross-validation step
                test_trials = random_set(x:(x + (ntrs_set - 1)));
                
                for d = 1:size(pooled_balanced_data,1)
                    
                    % reduce data set to the trials * number of sets (after dividing in X% sets)
                    temp_training_set = pooled_balanced_data{d,1,con}(:,:,1:(ntrs_set*cfg.cross_val_steps));
                    
                    if cfg.analysis_mode == 3 % if SVR
                        % reduce data set to the trials * number of sets (after dividing in X% sets)
                        % use specifiied varaible condition only
                        temp_training_labels = cfg.regress_data.SVR_matrix( 1:(ntrs_set*cfg.cross_val_steps) , cfg.regr_todo);
                    end % if
                    
                    % extract test-trials and delete them from training-trials
                    for trl = 1:size(test_trials,2)
                        
                        test_set{d,con,cv,ncv}(:,:,trl) = pooled_balanced_data{d,1,con}(:,:,test_trials(trl));
                        
                        if cfg.analysis_mode == 3 % if continuous SVR
                            % take the value for the same trial as data chosen for the test data set
                            cfg.test_labels{d,con,cv,ncv}(trl,1) = temp_training_labels(test_trials(trl),1);
                            cfg.test_trials{d,con,cv,ncv}(trl,1) = test_trials(trl);
                        end % if
                        
                        % option 1: delete used trials each for trial
                        % temp_training_set(:,:,delete_test_trials(trl))=[];
                        
                    end %trl
                    
                    % CROSS_CONDITION SVM__________________________________
                    % extract training set
                    if cfg.cross == 0
                        
                        delete_test_trials=fliplr(sort(test_trials));
                        
                    elseif cfg.cross == 1 || cfg.cross == 2
                        
                        % if doing cross-classification, ensure that trials
                        % from one training set don't end up in the
                        % opposite test set
                        delete_test_trials = [];
                        
                        for trl=1:size(test_trials,2)
                            
                            % find trials in training set which are the
                            % same as in opposite test set
                            for testTrl = 1:size(temp_training_set,3)
                                
                                if isequal(pooled_balanced_data{(abs(d-2)+1),1,con}(:,:,test_trials(trl)),temp_training_set(:,:,testTrl))
                                    
                                    delete_test_trials = [delete_test_trials testTrl];
                                    
                                end % if
                                
                            end % testTrl
                            
                        end %trl
                        
                        if numel(delete_test_trials) < numel(test_trials)
                            
                            delete_test_trials = [delete_test_trials, test_trials(1:(numel(test_trials) - numel(delete_test_trials)))];
                            
                        end % if
                        
                    end % cfg.cross if
                    %______________________________________________________
                    
                    % option 2: delete all used trials at once
                    temp_training_set(:,:,delete_test_trials) = [];
                    
                    training_set{d,con,cv,ncv} = temp_training_set;
                    
                    if cfg.analysis_mode == 3 % if continuous SVR
                        % delete the same trials from regression labels as for data
                        temp_training_labels(delete_test_trials) = [];
                        cfg.training_labels{d,con,cv,ncv} = temp_training_labels;
                    end
                    
                    clear delete_test_trials;
                    clear temp_training_set;
                    clear temp_training_labels;
                    
                    fprintf('Data sorted for single-trial analysis: specified DCG %d condition %d cross-validation step %d of cycle %d. \n',d,con,cv,ncv);
                    
                end % d
                
                % go to next step and repeat
                x=x+ntrs_set;
                
                clear test_trials
                
            end % cv (cross-validation steps)
            
            clear random_set
            
        end % ncv (repetition of cross-validation)
        
        clear ntrs_set;
        
    end % con (all conditions in serial order)
    
    clear pooled_balanced_data
    
% run-average data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% test- and training-data correspond to run averages 
elseif cfg.avmode == 2 % run averages
    
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

%% SECTION 7: BUILD LABELS AND DO CLASSIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

% input('Do you wish to continue with decoding? Press any button!')

fprintf('Starting with vector preparation... \n');
[RESULTS] = prepare_my_vectors_erp(training_set, test_set, cfg);

%% SECTION 8: AVERAGE CROSS-VALIDATION STEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

RESULTS.subj_acc = [];
RESULTS.subj_perm_acc = [];

% na = number analyses (stmode=2 has one per channel)
for na = 1:size(RESULTS.prediction_accuracy,2)
    
    pa(:,:,:) = RESULTS.prediction_accuracy{na}(:,:,:);
    
    if cfg.perm_test == 1
        perm_pa(:,:,:) = RESULTS.perm_prediction_accuracy{na}(:,:,:);
    end
               
    % calculate average decoding accuracy
    RESULTS.subj_acc(na,:) = nanmean(nanmean(pa,3),2);
    clear pa;
        
    % calculate average permutation test decoding accuracy
    if cfg.perm_test == 1
        RESULTS.subj_perm_acc(na,:) = nanmean(nanmean(perm_pa,3),2);
        clear perm_pa;
    end
        
end % na

fprintf('Results are computed and averaged for participant %d. \n',cfg.sbj);

%% SECTION 9: SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saves decoding results to a .mat file in the output directory. Some analysis
% settings (e.g. window_width_ms) are included in the file name.
%__________________________________________________________________________

if cfg.cross == 0
    savename = [(cfg.output_dir) cfg.study_name '_SBJ' num2str(cfg.sbj) '_win' num2str(cfg.window_width_ms) '_steps' num2str(cfg.step_width_ms)...
        '_av' num2str(cfg.avmode) '_st' num2str(cfg.stmode) '_' cfg.analysis_mode_label '_DCG' cfg.dcg_labels{cfg.dcg_todo} '.mat'];
elseif cfg.cross == 1
    savename = [(cfg.output_dir) cfg.study_name '_SBJ' num2str(cfg.sbj) '_win' num2str(cfg.window_width_ms) '_steps' num2str(cfg.step_width_ms)...
        '_av' num2str(cfg.avmode) '_st' num2str(cfg.stmode) '_' cfg.analysis_mode_label '_DCG' cfg.dcg_labels{cfg.dcg_todo(1)}...
        'toDCG' cfg.dcg_labels{cfg.dcg_todo(2)} '.mat'];
elseif cfg.cross == 2
        savename = [(cfg.output_dir) cfg.study_name '_SBJ' num2str(cfg.sbj) '_win' num2str(cfg.window_width_ms) '_steps' num2str(cfg.step_width_ms)...
        '_av' num2str(cfg.avmode) '_st' num2str(cfg.stmode) '_' cfg.analysis_mode_label '_DCG' cfg.dcg_labels{cfg.dcg_todo(2)}...
        'toDCG' cfg.dcg_labels{cfg.dcg_todo(1)} '.mat'];
end
    
save(savename,'cfg','RESULTS'); % Save cfg and RESULTS structures into a .mat file

fprintf('Results are saved for participant %d in directory: %s. \n',cfg.sbj,(cfg.output_dir));

%% SECTION 10: DISPLAY INDIVIDUAL RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displays decoding results for each individual subject dataset.
%__________________________________________________________________________

if cfg.display_on == 1
    
    display_indiv_results_erp(cfg, RESULTS);

end % display_on
%__________________________________________________________________________ 