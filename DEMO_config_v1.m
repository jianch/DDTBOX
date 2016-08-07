function DEMO_config_v1
%__________________________________________________________________________
% DDTBOX script written by Stefan Bode 01/03/2013
%
% The toolbox was written with contributions from:
% Daniel Bennett,  Daniel Feuerriegel, Phillip Alday
%
% The author further acknowledges helpful conceptual input/work from: 
% Jutta Stahl, Simon Lilburn, Philip L. Smith, Carsten Murawski, Carsten Bogler,
% John-Dylan Haynes
%__________________________________________________________________________
%
% This script is the configuration script for the DDTBOX. All
% study-specific information for decoding, regression and groupl-level
% analyses are specified here.
%
%__________________________________________________________________________
%
% Variable naming convention: STRUCTURE_NAME.example_variable

global SLIST;
global SBJTODO;

%% GENERAL STUDY PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

% Decide whether to save the SLIST structure and EEG data in a .mat file
savemode = 0;
ismymac = 2; % 1=mac; 2=PC
demo_mode = 2; % 1=results; 2=preanalysed; 3=SVR - for workshop demonstrations!

if ismymac==1

    bdir='/Volumes/PROJECTS/MVPA WORKSHOP COLOGNE/';
    
    if demo_mode == 1 
        output_dir='/Volumes/PROJECTS/MVPA WORKSHOP COLOGNE/DECODING_RESULTS/results/';
    elseif demo_mode == 2
        output_dir='/Volumes/PROJECTS/MVPA WORKSHOP COLOGNE/DECODING_RESULTS/preanalysed/';
    elseif demo_mode == 3
        output_dir='/Volumes/PROJECTS/MVPA WORKSHOP COLOGNE/DECODING_RESULTS/SVR/';
    end
    
    % subject codes/names
    sbj_code = {...

        ['DATA/sbj1/SBJ1_full'];... %1
        ['DATA/sbj2/SBJ2_full'];... %2 
        ['DATA/sbj3/SBJ3_full'];... %3
        ['DATA/sbj4/SBJ4_full'];... %4
        ['DATA/sbj5/SBJ5_full'];... %5

        };
    
elseif ismymac==2

    bdir='F:\MVPA_WORKSHOP\';
    
    if demo_mode == 1 
        output_dir='F:\MVPA_WORKSHOP\DECODING_RESULTS\results\';
    elseif demo_mode == 2    
        output_dir='F:\MVPA_WORKSHOP\DECODING_RESULTS\preanalysed\';
    elseif demo_mode == 3    
        output_dir='F:\MVPA_WORKSHOP\DECODING_RESULTS\SVR\';
    end

    % subject codes/names
    sbj_code = {...

        ['DATA\sbj1\SBJ1_full'];... %1
        ['DATA\sbj2\SBJ2_full'];... %2 
        ['DATA\sbj3\SBJ3_full'];... %3
        ['DATA\sbj4\SBJ4_full'];... %4
        ['DATA\sbj5\SBJ5_full'];... %5

        };   
    
end

nsbj=size(sbj_code,1);


%% CREATE SLIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

SLIST = [];
sn = SBJTODO;
   
    % subject parameters
    SLIST.number = sn;
    SLIST.sbj_code = sbj_code{sn};
    SLIST.output_dir = output_dir;
    SLIST.data_struct_name = 'eeg_sorted_cond';
    
    % if SVR
    SLIST.regress_label_name = [bdir sbj_code{sn} 'regress_sorted_data.mat'];
    SLIST.regress_struct_name='SVR_matrix'; % DO NOT CHANGE NAME
    
    % channels    
    SLIST.nchannels=64;
    SLIST.channels='channel_labels';
    SLIST.channel_names_file='channel_inf.mat';
    SLIST.channellocs=[bdir 'locations/'];
    
    % sampling rate and baseline
    SLIST.sampling_rate=1000;
    SLIST.pointzero=100; % corresponds to zero (time-locked to this event, in ms)
        
%% CREATE DCGs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    
%% SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

% Save the SLIST structure and eeg_sorted_cond to a .mat file
if savemode == 1
    
    % DF NOTE: I have changed the second argument from 'eeg_sorted_cond' to
    % SLIST.data_struct_name so that it will still save the EEG data file
    % if the user decides to use a different variable name than
    % 'eeg_sorted_cond'
    save(SLIST.data_save_name, SLIST.data_struct_name, 'SLIST');
    
end  

