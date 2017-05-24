function dd_convert_eeg_data(EEG, events_by_cond, save_directory, save_filename, eeg_toolbox, data_type)
% dd_convert_epoched_eeg_data.m
%
% This function extracts EEG data epoched in EEGLAB or
% ERPLAB for use in DDTBOX, and saves the DDTBOX-compatible EEG data in the
% following array:
%
%   eeg_sorted_cond{run, condition}(timept, channel, epoch)
%
% Epoching and artefact rejection must be completed prior to running this
% function, either in EEGLAB or ERPLAB. If using EEGLAB, make sure to
% remove epochs containing artefacts from your data first, by going to
% Tools -> Reject data epochs -> Reject marked epochs
%
% WARNING: This script is a beta version and has not been thoroughly
% tested across versions of EEGLab/ERPLab and on different EEG data types.
% Wherever possible, always check whether this script is properly
% extracting the correct epochs by manually extracting epoched data from 
% the EEG structure in EEGLab. The code in this function should give you an
% idea of how to do this.
%
%
% Inputs:
%   
%   EEG     Structure containing EEG data and epoch information, created
%           EEGLAB.
%
%   events_by_cond      Array containing event codes (if using EEGLAB) or bin
%                       indices (if using ERPLAB) for epochs in each
%                       condition for MVPA.
%
%   save_directory       Filepath for saving the resulting DDTBOX-compatible
%                       .mat file containing epoched EEG data.
%
%   save_filename       Name of the .mat file containing DDTBOX-compatible
%                       epoched EEG data.
%
%   eeg_toolbox         Name of the toolbox used for epoching EEG data.
%                       This function accepts either 'EEGLAB' or 'ERPLAB'.
% 
%   data_type           Select whether to extract EEG data or independent component activations.
%                       This function accepts either 'EEG' or 'ICAACT'
%
% Optional keyword inputs:
%
% 
% Example:  dd_convert_eeg_data(EEG, 'DDTBOX-Data/ID1/', 'ID1.mat', 'EEGLAB', 'ICAACT')
%
% 
%
%
% Copyright (c) 2013-2017, Daniel Feuerriegel and contributors 
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



% Determine data type to use
if strcmp(data_type, 'EEG') == 1 % If using EEG data
    epoched_data = EEG.data;
    n_epochs_total = size(epoched_data, 3);
elseif strcmp(data_type, 'ICAACT') == 1 % If using independent component activations
    epoched_data = EEG.icaact;
    n_epochs_total = size(epoched_data, 3);
else % If data type not correctly specified
    error('Data type for DDTBOX not correctly specified. Please input either "EEG" or "ICAACT"');
end % of if strcmp data_type


% Check if save filepath exists, create directory if doesn't exist
if ~exist(save_directory, 'dir')
    fprintf(['\n\n Specified save directory does not exist.\n Creating the directory' save_directory '...\n\n']);
    mkdir(save_directory);
end % of if ~exist

% Check number of conditions and number of event codes in each condition
n_conds = length(events_by_cond);

n_eventcodes_by_cond = nan(1, n_conds); % Preallocate
for cond_no = 1:n_conds
    n_eventcodes_by_cond(cond_no) = length(events_by_cond{cond_no});
end % of for cond_no

% Create empty eeg_sorted_cond cell array
eeg_sorted_cond = cell(1, n_conds);


% Different methods of extracting epochs are used depending on whether data
% was epoched in EEGLAB or ERPLAB
if strcmp(eeg_toolbox, 'EEGLAB') == 1 % If using EEGLAB
    
    for epoch_no = 1:n_epochs_total
        
        bin_indices = (EEG.epoch(epoch_no).eventtype);
    
        % Checking whether event codes are strings and converting to double
        % precision floating point
        if ischar(bin_indices{1})
            for bin_index_temp = 1:length(bin_indices)
                bin_indices{bin_index_temp} = str2num(bin_indices{bin_index_temp});
            end % of for bin_index_temp
        end % of if isstring
        
        % Convert cell array to vector
        bin_indices = cell2mat(bin_indices);
        
        for bin_index = 1:length(bin_indices) % For each bin index in the epoch
        
            % Cycle through each condition and look for matching event codes
            for condition_no = 1:n_conds
                for event_code_no = 1:n_eventcodes_by_cond(condition_no)
                    
                    % Check whether bin index matches a specified event
                    % code corresponding to a condition in DDTBOX analyses
                    if bin_index == events_by_cond{condition_no}(event_code_no)
                        
                         % Check if artefact rejection has been conducted, and
                         % if the epoch has been marked for rejection
                         % using artefact detection routines
                         if ~isempty(EEG.reject.rejmanual); % If artefact detection has been conducted
                             if EEG.reject.rejmanual(1, epoch_no) == 0 % If not marked for rejection

                                 % Copy the epoch into the eeg_sorted_cond cell array
                                 eeg_sorted_cond{1, condition_no}(:,:,end + 1) = epoched_data(:,:,epoch_no);

                             end % of if EEG.reject.rejmanual
                         end % of if ~isempty EEG.reject.rejmanual
                    end % of if bin_index 
                end % of for event_code_no
            end % of for condition_no 
        end % of for bin_index
    
    
    % Note: EEGLAB uses EEG.reject.rejmanual like ERPLAB so we can do
    % automatic checking like in ERPLAB to automatically detect and exclude
    % bad epochs.
    
    
    end % of for epoch_no
    
    
elseif strcmp(eeg_toolbox, 'ERPLAB') == 1 % If using ERPLab
    
    % Go through each epoch and check whether it belongs to a bin index
    % specified in events_by_cond
    for epoch_no = 1:n_epochs_total % Go through all epochs
        
        % Get vectors of bin indices for the epoch
        bin_indices = cell2mat(EEG.epoch(epoch_no).eventbini);
        
        for bin_index = 1:length(bin_indices) % For each bin index in the epoch
        
            % Cycle through each condition and look for matching event codes
            for condition_no = 1:n_conds
                for event_code_no = 1:n_eventcodes_by_cond(condition_no)
                    
                    % Check whether bin index matches a specified event
                    % code corresponding to a condition in DDTBOX analyses
                    if bin_index == events_by_cond{condition_no}(event_code_no)
                        
                         % Check if artefact rejection has been conducted, and
                         % if the epoch has been marked for rejection
                         % using artefact detection routines
                         if ~isempty(EEG.reject.rejmanual); % If artefact detection has been conducted
                             if EEG.reject.rejmanual(1, epoch_no) == 0 % If not marked for rejection

                                 % Copy the epoch into the eeg_sorted_cond cell array
                                 eeg_sorted_cond{1, condition_no}(:,:,end + 1) = epoched_data(:,:,epoch_no);

                             end % of if EEG.reject.rejmanual
                         end % of if ~isempty EEG.reject.rejmanual
                    end % of if bin_index 
                end % of for event_code_no
            end % of for condition_no 
        end % of for bin_index
    end % of for epoch_no
    
else % EEG toolbox name incorrectly specified
    
    error('EEG toolbox name not correctly specified. Please input either "EEGLAB" or "ERPLAB"');

end % of if strcmp eeg_toolbox


% Save resulting DDTBOX-compatible epoched EEG data file
try
    save([save_directory, save_filename], 'eeg_sorted_cond');
catch % If user has forgotten forward slash
    save([save_directory, '/', save_filename], 'eeg_sorted_cond');
end


