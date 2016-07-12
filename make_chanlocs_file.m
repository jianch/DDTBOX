function make_chanlocs_file(EEG, varargin)

%__________________________________________________________________________
% Script to make DDTBox channel location files from a loaded EEGLab dataset. Written by
% Daniel Feuerriegel on 11/07/2016.
%
% The toolbox was written with contributions from:
% Daniel Bennett, Daniel Feuerriegel, Phillip Alday
%
% The author (Stefan Bode) further acknowledges helpful conceptual input/work from: 
% Jutta Stahl, Simon Lilburn, Philip L. Smith, Elaine Corbett, Carsten Murawski, 
% Carsten Bogler, John-Dylan Haynes
%__________________________________________________________________________
%
% This script copies the channel location information from a loaded EEGLab
% dataset and uses this data to create a channel locations file. The channel locations file
% is then saved at the specified location. An EEGLab dataset must be loaded 
% first for this script to work.
%
% required:
% - EEG (EEGLab data structure)
% 
% optional:
% - save_filepath (location in which to save the chanlocs file) 
%
%__________________________________________________________________________
%
% Variable naming convention: STRUCTURE_NAME.example_variable

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
