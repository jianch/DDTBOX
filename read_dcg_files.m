function [DCG] = read_dcg_files(file_path)
% This function reads text files that describe discrimination groups used
% in DDTBOX and outputs a structure containing the discrimination group
% name, condition names and event codes for epochs associated with each
% condition. 
%
%
% The organisation of dcg files is as follows:
% DCG Name
% Condition 1 name
% Condition 1 event codes
% Condition 2 name
% Condition 2 event codes
%
% Example:
% Faces vs. Objects
% Faces
% 10, 11, 12, 13
% Objects
% 20, 21, 22, 23
%
%
%
% Inputs:
%
%   file_path   The file path of the text file containing discrimination
%               group information.
%
%
% Outputs:
%
%   DCG Structure containing:
%
%   dcg_name            The name of the discrimination group
%
%   cond1_name          The label for condition 1
%
%   cond1_event_codes   Vector of event codes of epochs associated with condition 1
%
%   cond2_name          The label for condition 2
%
%   cond2_event_codes   Vector of event codes of epochs associated with condition 2
%
%
% Example:      [DCG] = read_dcg_files('dcg files/current study/example_dcg.txt')
%
%
%
% Copyright (c) 2016 Daniel Feuerriegel and contributors
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


% Open the text file containing dcg information
file_ID = fopen(file_path);

% Copy the dcg text file into a cell array
temp_dcg_info = textscan(file_ID,'%s', 5, 'Delimiter', '\n');

% Copy information into the DCG structure

% DCG and condition names
DCG.dcg_name = char(temp_dcg_info{1, 1}(1)); % Discrimination group name
DCG.cond1_name = char(temp_dcg_info{1, 1}(2)); % Condition 1 name
DCG.cond2_name = char(temp_dcg_info{1, 1}(4)); % Condition 4 name

% Event codes for each condition

% Condition 1
temp_string = char(temp_dcg_info{1, 1}(3)); % Convert cell to string
temp_event_codes = textscan(temp_string ,'%u16', 'Delimiter', ','); % Condition 1
DCG.cond1_event_codes = cell2mat(temp_event_codes);

% Condition 2
temp_string = char(temp_dcg_info{1, 1}(5)); % Convert cell to string
temp_event_codes = textscan(temp_string ,'%u16', 'Delimiter', ','); % Condition 1
DCG.cond2_event_codes = cell2mat(temp_event_codes);


end % of function
