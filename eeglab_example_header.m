function eeglab_example_header(EEG, varargin)
%
% This script copies the channel location information from a loaded EEGLab
% dataset and uses this data to create a channel locations file. The channel locations file
% is then saved at the specified location. An EEGLab dataset must be loaded 
% first for this script to work.
%
%
% Required Inputs:
% - EEG (EEGLab data structure)
% 
% Optional Inputs:
% - save_filepath (location in which to save the chanlocs file) 
%
% Outputs:
% 
% Example:
%
%
% Copyright (c) 2016 _________ and contributors
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
