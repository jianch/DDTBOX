function EXAMPLE_config_v1
%
% This is a configuration script for the DDTBOX. All
% study-specific information for decoding, regression and group-level
% analyses are specified here. 
%
% Copyright (c) 2013-2016 Stefan Bode and contributors
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


% Set SLIST and SBJTODO to globals
global SLIST;
global SBJTODO;

%% GENERAL STUDY PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Configurable settings are:
%   savemode            Decide whether to save a copy of the configuration
%                       parameters in a .mat file
%
%   bdir                Base directory path (where single subject EEG datasets
%                       will be stored)
% 
%   output_dir          Output directory (where decoding results will be
%                       saved)
%
%   sbj_code            Filepaths of single subject datasets (relative to 
%                       the base directory path)

%__________________________________________________________________________

% Decide whether to save the SLIST structure .mat file (to keep a separate
% record of your configuration settings)
savemode = 0; % 1 = save / 0 = don't save

% Base directory path (where single subject EEG datasets will be stored)
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
    

% Calculates number of subjects from the number of data files
nsbj = size(sbj_code, 1);


%% CREATE SLIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The SLIST structure contains all study parameters set in the
% configuration file. In this section we set EEG recording parameters and
% the names of the MATLAB single subject data arrays.
%
% Configurable settings are:
%   SLIST.data_struct_name      The name of single subject data arrays in
%                               the MATLAB workspace.
%
%   SLIST.regress_label_name    The filepath of .mat files containing condition 
%                               labels for support vector regression.
%
%   SLIST.nchannels             The number of channels in the EEG dataset
%
%   SLIST.channels              ???
%
%   SLIST.channel_names_file    The filename of the .mat file containing
%                               EEG channel locations
%
%   SLIST.channellocs           The path of the folder of the EEG channel 
%                               locations .mat file.
%
%   SLIST.sampling_rate         The sample rate of the EEG data (in Hz)
%
%   SLIST.pointzero             The time of the event/trigger code relative
%                               to the beginning of the epoch (in ms). For
%                               a stimulus-locked response this would be
%                               the length of the prestimulus baseline
%                               period.
%
%__________________________________________________________________________

SLIST = [];
sn = SBJTODO;
   
    % Subject parameters
    SLIST.number = sn;
    SLIST.sbj_code = sbj_code{sn};
    SLIST.output_dir = output_dir;
    SLIST.data_struct_name = 'eeg_sorted_cond'; % Data arrays for use with DDTBOX must use this name as their MATLAB workspace
    
    % Settings for support vector regression
    SLIST.regress_label_name = [bdir, sbj_code{sn}, 'regress_sorted_data.mat'];
    SLIST.regress_struct_name = 'SVR_matrix'; % DO NOT CHANGE NAME
    
    % channels    
    SLIST.nchannels = 64;
    SLIST.channels = 'channel_labels';
    SLIST.channel_names_file = 'channel_inf.mat';
    SLIST.channellocs = [bdir, 'locations/'];
    
    % sampling rate and baseline
    SLIST.sampling_rate = 1000; % Data sample rate in Hz
    SLIST.pointzero = 100; % Corresponds to the time of the event relative to the prestimulus baseline (in ms)
        
    
%% CREATE DCGs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
       
    SLIST.ndcg = size(SLIST.dcg,2);
    SLIST.nclasses = size(SLIST.dcg{1},2);      
    SLIST.ncond = size(SLIST.cond_labels,2);
    
    SLIST.data_open_name = [bdir (sbj_code{sn}) '.mat'];
    SLIST.data_save_name = [bdir (sbj_code{sn}) '_data.mat'];
    
    
%% SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

% Save the SLIST structure and eeg_sorted_cond to a .mat file
if savemode == 1
    
    save(SLIST.data_save_name, SLIST.data_struct_name, 'SLIST');
    
end % of if savemode

