% decoding_erp() - Perform MVPA on a single subject 
%
% Usage: 
%		 >> [results, cfg] = decoding_erp(dat,cfg,'Key1',Value1,...);
%
% Inputs:
%
%   dat            Either a file name (including path) of file containing
%                  floating point data, or a data matrix (chans x frames)
%   cfg			   A MATLAB struct containing the necessary configuration, 
%				   see pop_dd_decode() for more information on necessary fields				   
%  'Key1'          Keyword string for argument 1
%   Value1         Value of argument 1
%   ...            ...
%		
% Optional keyword inputs:
%
%   outdir         name of directory to write output (does not have to exist), def=pwd/ddouttmp/
%   indir          optional input directory from which to load init
%
% Outputs:
%
%   results        classification results (for use in e.g. analyse_decoding_erp())
%   cfg            updated configuration with 'dependent' parameters calculated from 
% 				   'independent' parameters, e.g window_width (in samples) calculated from 
%				   window_width_ms (in milliseconds)
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

function [results, cfg] = decoding_erp(dat, cfg, varargin)
cfg.window_width = floor(cfg.window_width_ms / ((1/cfg.srate) * 1000));
cfg.step_width = floor(cfg.step_width_ms / ((1/cfg.srate) * 1000));

% cross-validation defaults for single-trial decoding
cfg.n_all_analyses = cfg.cross_val_steps * cfg.n_rep_cross_val;
cfg.n_all_permutation = cfg.cross_val_steps * cfg.permut_rep;

if cfg.cross == 0 
    fprintf('DCG %d will be analysed. \n',cfg.dcg_todo);                
elseif cfg.cross > 0
    fprintf('DCG %d and %d will be analysed for cross-condition classification. \n',cfg.dcg_todo(1), cfg.dcg_todo(2));
end
    
fprintf('Cross-validation defaults for single-trial decoding:\n');
fprintf('%d steps with %d cycles resulting in %d analyses.\n',cfg.cross_val_steps,cfg.n_rep_cross_val,cfg.n_all_analyses);
fprintf('Balanced number of examples will be used for all analyses by default.\n');

if cfg.perm_test == 1
    fprintf('Random-label analysis will be based on %d analyses.\n',cfg.n_all_permutation);
end

%% SECTION 2: READ IN DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic data in: eeg_sorted_cond{run,cond}(timepoints,channels,trials)
% *** converted into work_data{run,cond}(timepoints,channels,trials)

fprintf('Reading in data. Please wait... \n');
open_name = (SLIST.data_open_name);
load(open_name);
fprintf('Data loading complete.\n');

% read in regression labels if performing continuous regression
if cfg.analysis_mode == 4 %i.e. if we are performing continuous regression
    
   cfg.training_label_open_name = SLIST.training_label_file;
   cfg.training_label_import = load(cfg.training_label_open_name);
   cfg.test_label_open_name = SLIST.test_label_file;
   cfg.test_label_import = load(cfg.test_label_open_name);
   fprintf('Loaded regressand labels.\n');

end

%__________________________________________________________________________
% if rt_match==1
%     
%     % match_my_rt.m is a toolbox-script that finds RT-matched trials in two
%     % conditions form dcg_todo. Only these conditions are matched, the
%     % other condition's data is carried forward but not used for decoding.
%     fprintf('Peforming RT-match for conditions.\n');
%     [eeg_sorted_cond_matched_rt,ntrials_dcg_todo,final_rts]=match_my_rts(eeg_sorted_cond,rt_sorted_cond,cfg.dcg_todo);
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
% cfg_slist defines all DCG with their respective relevant conditions. In
% this section, work_data is reduced to the specified conditions
% *** work_data will be reduced to reduced_data{dcg,run,cond}(timepoints,channels,trials)
% (this section was modified in v29 to accomodate cross-condition decoding. A dimension was added to the data-matrix
%__________________________________________________________________________

if cfg.analysis_mode ~= 4 % continuous SVR does not require conditions
          
    for r = 1:size(work_data,1) % for all runs
        
        % go through either one DCG (regular) or two DCGs (cross-condition decoding)
        for d = 1:size(cfg.dcg_todo,2)
            
            for cond = 1:size(SLIST.dcg{cfg.dcg_todo(d)},2) % for all conditions specified

                fprintf('Run %d: Extracting condition %d as specified in DCG %d.\n',r,(SLIST.dcg{cfg.dcg_todo(d)}(cond)),cfg.dcg_todo(d));
                temp(:,:,:) = work_data{r,(SLIST.dcg{cfg.dcg_todo(d)}(cond))}(:,:,:);
                reduced_data{d,r,cond} = temp; 

                clear temp;                    

            end % cond    
        end % d
    end % r 

elseif cfg.analysis_mode == 4 % put in by Dan 14/3/2014
    reduced_data{1,:,:} = work_data{:,:}; 
end % of if cfg.analysis_mode ~= 4 statement

cfg.nconds = size(reduced_data,3); % Calculates the number of conditions
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
            cfg.nextra = SLIST.extra;

            if cfg.nextra > 0
                fprintf('Removing extra channels for participant %d condition %d \n',SBJTODO,cond);
                extra_channels = SLIST.extra;
                temp_data(:,extra_channels,:) = [];
            elseif cfg.nextra == 0      
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
    cfg.mintrs_min(r) = min(mintrs);
    
end %r

fprintf('Minimum number of trials per condition computed for participant %d \n',SBJTODO);

% select min trials in all conditions
for r = 1:size(clean_data,2)
    
    for d = 1:size(clean_data,1)

        for cond = 1:size(clean_data,3)

            % extract min number trials for this DCG and run
            temp_data = clean_data{d,r,cond};
            balanced_data{d,r,cond} = temp_data(:,:,1:(cfg.mintrs_min(r)));
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


%% SECTION 7: SORT DATA FOR CLASSIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                    
                    temp_training_set = pooled_balanced_data{d,1,con}(:,:,:);
                    
                    if cfg.analysis_mode == 4 % if continuous SVR
                        temp_training_labels = cfg.training_label_import.labels{1,con};
                    end % if
                    
                    % extract test-trials and delete them from training-trials
                    for trl = 1:size(test_trials,2)
                        
                        test_set{d,con,cv,ncv}(:,:,trl) = pooled_balanced_data{d,1,con}(:,:,test_trials(trl));
                        
                        if cfg.analysis_mode == 4 % if continuous SVR
                            cfg.test_labels{con,cv,ncv}(trl) = cfg.test_label_import.labels{1,con}(test_trials(trl));
                        end % if
                        
                        % option 1: delete used trials each for trial
                        % temp_training_set(:,:,delete_test_trials(trl))=[];
                        
                    end %trl
                    
                    
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
                    
                    % option 2: delete all used trials at once
                    temp_training_set(:,:,delete_test_trials) = [];
                    
                    training_set{d,con,cv,ncv} = temp_training_set;
                    
                    if cfg.analysis_mode == 4 % if continuous SVR
                        temp_training_labels(delete_test_trials) = [];
                        cfg.training_labels{con,cv,ncv} = temp_training_labels;
                    end
                    
                    clear delete_test_trials;
                    clear temp_training_set;
                    clear temp_training_labels;
                    
                    fprintf('Data sorted for single-trial decoding: specified DCG %d condition %d cross-validation step %d of cycle %d. \n',d,con,cv,ncv);
                    
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

%% SECTION 8: BUILD LABELS AND DO CLASSIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

% input('Do you wish to continue with decoding? Press any button!')

fprintf('Starting with vector preparation... \n');
[RESULTS] = prepare_my_vectors_erp(training_set, test_set, SLIST, cfg);

%% SECTION 9: AVERAGE CROSS-VALIDATION STEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% SECTION 10: SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saves decoding results to a .mat file in the output directory. Some analysis
% settings (e.g. window_width_ms) are included in the file name.
%__________________________________________________________________________

if cfg.cross == 0
    savename = [(SLIST.output_dir) cfg.cfg_name '_SBJ' num2str(cfg.sbj) '_win' num2str(cfg.window_width_ms) '_steps' num2str(cfg.step_width_ms)...
        '_av' num2str(cfg.avmode) '_st' num2str(cfg.stmode) '_DCG' SLIST.dcg_labels{cfg.dcg_todo} '.mat'];
elseif cfg.cross == 1
    savename = [(SLIST.output_dir) cfg.cfg_name '_SBJ' num2str(cfg.sbj) '_win' num2str(cfg.window_width_ms) '_steps' num2str(cfg.step_width_ms)...
        '_av' num2str(cfg.avmode) '_st' num2str(cfg.stmode) '_DCG' SLIST.dcg_labels{cfg.dcg_todo(1)} 'toDCG' SLIST.dcg_labels{cfg.dcg_todo(2)} '.mat'];
elseif cfg.cross == 2
        savename = [(SLIST.output_dir) cfg.cfg_name '_SBJ' num2str(cfg.sbj) '_win' num2str(cfg.window_width_ms) '_steps' num2str(cfg.step_width_ms)...
        '_av' num2str(cfg.avmode) '_st' num2str(cfg.stmode) '_DCG' SLIST.dcg_labels{cfg.dcg_todo(2)} 'toDCG' SLIST.dcg_labels{cfg.dcg_todo(1)} '.mat'];
end
    
save(savename,'cfg','RESULTS'); % Save cfg and RESULTS structures into a .mat file

fprintf('Results are saved for participant %d in directory: %s. \n',cfg.sbj,(SLIST.output_dir));

%% SECTION 11: DISPLAY INDIVIDUAL RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displays decoding results for each individual subject dataset.
%__________________________________________________________________________

if cfg.display_on == 1
    
    display_indiv_results_erp(cfg,RESULTS);

end % display_on
%__________________________________________________________________________ 