% EXAMPLE_run_decoding_analyses.m
%
% This script is used for specifying configuration settings for DDTBOX and
% running decoding analyses on specified subjects and discrimination
% groups. An explanation of each configurable parameter is described below.
% Please make copies of this script for your own projects.
%
% This script calls decoding_erp.m
%
%
% Copyright (c) 2013-2016, Daniel Feuerriegel and contributors 
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


%% Filepaths and locations of subject datasets

% Base directory path (where single subject EEG datasets and channel locations information will be stored)
bdir = '/Volumes/PROJECTS/MVPA WORKSHOP COLOGNE/';

% Output directory (where decoding results will be saved)
output_dir = '/Volumes/PROJECTS/MVPA WORKSHOP COLOGNE/DECODING_RESULTS/results/';
    
% Filepaths of single subject datasets (relative to the base directory)
sbj_code = {...

    ['DATA/sbj1/SBJ1_full'];... %subject 1
    ['DATA/sbj2/SBJ2_full'];... %subject 2 
    ['DATA/sbj3/SBJ3_full'];... %subject 3
    ['DATA/sbj4/SBJ4_full'];... %subject 4
    ['DATA/sbj5/SBJ5_full'];... %subject 5

    };
    

% Automatically calculates number of subjects from the number of data files
nsbj = size(sbj_code, 1);


% MATLAB workspace names for single subject data arrays and structures
SLIST.data_struct_name = 'eeg_sorted_cond'; % Data arrays for use with DDTBOX must use this name as their MATLAB workspace
    
% Settings for support vector regression
SLIST.regress_label_name = [bdir, sbj_code{sn}, 'regress_sorted_data.mat'];
SLIST.regress_struct_name = 'SVR_matrix'; % DO NOT CHANGE NAME



%% EEG dataset information

nchannels = 64; % number of channels
channel_names_file = 'channel_inf.mat'; % Name of .mat file containing channel labels and channel locations
channellocs = [bdir, 'locations/']; % Path of directory containing channel information file
sampling_rate = 1000; % Data sampling rate in Hz
pointzero = 100; % Corresponds to the time of the event/trigger code relative to the prestimulus baseline (in ms)



%% Condition and discrimination group information
%
% Discrimination Groups (DCGs) are used to label and organise conditions to
% be discriminated between using MVPA.
%
% Configurable settings are:
%   SLIST.cond_labels{X}        Names of each condition. Conditon label {X}
%                               corresponds to data in column X of the
%                               single subject data arrays.
%
%   SLIST.dcg                   Indices of conditions to be discriminated
%                               between for each DCG. For example, to
%                               discriminate between conditions 3 and 4 for
%                               DCG 1 we input: SLIST.dcg{1} = [3, 4];  
%                               
%   SLIST.dcg_labels            The name/label for each DCG. For example
%                               'Faces vs. Chairs' or 
%                               'Errors vs. Correct Responses'
%
%__________________________________________________________________________

% Label each condition
% Example: SLIST.cond_labels{condition number} = 'Name of condition';
SLIST.cond_labels{1} = 'condition_A';
SLIST.cond_labels{2} = 'condition_B';
SLIST.cond_labels{3} = 'condition_C';
SLIST.cond_labels{4} = 'condition_D';
        
% Discrimination groups
% Enter the condition numbers of the conditions to discriminate between
% Example: SLIST.dcg{discrimination group number} = [condition number 1, condition number 2];
SLIST.dcg{1} = [1, 3]; 
SLIST.dcg{2} = [2, 4]; 
              
% Label each discrimination group
% Example: SLIST.dcg_labels{Discrimination group number} = 'Name of discrimination group'
SLIST.dcg_labels{1} = 'A vs. B';
SLIST.dcg_labels{2} = 'C vs. D';



% This section automaticallly fills in various parameters related to dcgs and conditions 
SLIST.ndcg = size(SLIST.dcg, 2);
SLIST.nclasses = size(SLIST.dcg{1}, 2);      
SLIST.ncond = size(SLIST.cond_labels, 2);

SLIST.data_open_name = [bdir, (sbj_code{sn}), '.mat'];
SLIST.data_save_name = [bdir, (sbj_code{sn}), '_data.mat'];








%% Select experiment, subject datasets and discrimination groups

% Enter the name of the study (used to find the study config script file). 
% e.g. if the study_name = 'EXAMPLE' and vconf = 1 then DDTBOX will use
% the config script EXAMPLE_config_v1.m
study_name = 'EXAMPLE';
vconf = 1;

% Set which subjects datasets to decode
sbj = [];

% Enter the discrimination group for classification. Two discrimination
% groups can be entered when using cross-condition decoding.
dcg_todo = [];

% Perform cross-condition decoding? 1 = yes / 0 = no
cross = 0;


%% Multivariate classification/regression parameters

analysis_mode = 1; % ANALYSIS mode (1=SVM with LIBSVM / 2=SVM with liblinear / 3=SVR with LIBSVM)
stmode = 1; % SPACETIME mode (1=spatial / 2=temporal / 3=spatio-temporal)
avmode = 1; % AVERAGE mode (1=no averaging; single-trial / 2=run average) !!! need single trials for SVR !!!
window_width_ms = 10; % width of sliding window in ms
step_width_ms = 10; % step size with which sliding window is moved through the trial
zscore_convert = 0; % Convert data into z-scores before decoding? 0 = no / 1 = yes
perm_test = 1; % Run decoding using permuted condition labels? 0=no / 1=yes
cross_val_steps = 5; % How many cross-validation steps (if no runs available)?
n_rep_cross_val = 10; % How many repetitions of full cross-validation with re-ordered data?
permut_rep = 10; % How many repetitions of full cross-validation with permutation results?

% Feature weights extraction
feat_weights_mode = 1; % Extract feature weights? 0=no / 1=yes

% Single subject decoding results plotting
display_on = 1; % Display individual subject results? 1=figure displayed / 0=no figure
perm_disp = 1; % display the permutation decoding results in figure? 0=no / 1=yes



%% Run the decoding analyses

for dcg = 1:length(discriminationGroups)
    for sbj = subjects

        cfg.sbj = sbj;
        cfg.dcg = dcg;
        
        DECODING_ERP(cfg);

    end % of for sbj
end % of for dcg