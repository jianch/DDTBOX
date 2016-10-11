function make_chanlocs_file(EEG, varargin)
%
% This script copies the channel location information from a loaded EEGLab
% dataset and uses this data to create a channel locations file. The channel locations file
% is then saved at the specified location. An EEGLab dataset must be loaded 
% first for this script to work.
%
% Inputs:
%
%   EEG             EEGLab data structure
% 
% optional:
%   save_filepath 	location in which to save the channel locations file 
%
% Outputs:
%
% Example:          make_chanlocs_file(EEG, 'Channel Locations/chanlocs_for_DDTBox.mat');
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

%% Handling variadic inputs
% Define defaults at the beginning
options = struct(...
    'save_filepath', []);

% Read the acceptable names
option_names = fieldnames(options);

% Count arguments
n_args = length(varargin);
if round(n_args/2) ~= n_args/2
   error([mfilename ' needs property name/property value pairs'])
end

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
   inp_name = lower(pair{1}); % make case insensitive

   % Overwrite default options
   if any(strcmp(inp_name, option_names))
      options.(inp_name) = pair{2};
   else
      error('%s is not a recognized parameter name', inp_name)
   end
end
clear pair
clear inp_name

% Renaming variables for use below:
save_filepath = options.save_filepath;
clear options;


%% Creating the chanlocs file

% Copy the relevant channel information from the EEGLab data structure
chaninfo = EEG.chaninfo;
chanlocs = EEG.chanlocs;

% Saves the file in the specified location
save([save_filepath 'chanlocs.mat'], 'chaninfo', 'chanlocs');
