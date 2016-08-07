% pop_dd_decode() - Convenience interface for calling decoding_erp()
%
% Usage: 
%		 >> [EEG, com] = pop_dd_decode(EEG,cfg,'Key1',Value1,...);
%
% Inputs:
%
%   dat            Either a file name (including path) of file containing
%                  floating point data, or a data matrix (chans x frames)	   
%  'Key'          Keyword string 
%   Value         Value
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
% Copyright (c) 2016, Phillip Alday and contributors 
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

function [EEG, com] = pop_dd_decode(EEG,varargin)

com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            

          
% defaults to populate the GUI with
% note that these defaults can't be used without popping up the GUI -- this
% is intentional as these may change. If you like them, then use the eegh
% mechanism to copy the final call generated with them all filled in.

cfg.stmode = DDSTMode.spatial; % SPACETIME mode 
cfg.avmode = DDAverage.none; % AVERAGE mode 
cfg.analysis_mode = DDMethod.svm; % ANALYSIS mode 
cfg.backend = DDBackend.libsvm;

cfg.window_width_ms = 10; % width of sliding window in ms
cfg.step_width_ms = 10; % step size with which sliding window is moved through the trial

cfg.perm_test = 1; % run the permutation-decoding? 
cfg.perm_disp = 1; % display the permutation results in figure? 
cfg.display_on = 1; % display figure?
cfg.pointzero = 'auto'; % where does the zero point fall within the epoch? 
% for e.g. an epoch -200..1200, the time locking eventing is 200ms from epoch start

cfg.rt_match = 0; % Use RT-matching algorithm to select trials? 
cfg.zscore_convert = 0; % Convert data into z-scores before decoding? 
cfg.feat_weights_mode = 1; % Extract feature weights? 
cfg.cross_val_steps = 10; % How many cross-validation steps (if no runs available)?
cfg.n_rep_cross_val = 10; % How many repetitions of full cross-validation with randomly re-ordered data?
cfg.permut_rep = 10; % How many repetitions of full cross-validation with permutation results?

% These flags are not currently exposed in the GUI or pop_() interface 
% because it's assumed that if you're advanced enough to be using these, 
% you're advanced enough to be calling these methods directly
backend_flags.svm_type = 0;
backend_flags.kernel_type = 0;
backend_flags.cost = 1;
backend_flags.extra_flags = [];

% Read the acceptable names
optionNames = fieldnames(cfg);

% DCGs


% Count arguments
nArgs = length(varargin);
if round(nArgs/2) ~= nArgs/2
   error([mfilename ' needs property name/property value pairs'])
end

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
   inpName = lower(pair{1}); % make case insensitive

   % Overwrite default options
   if any(strcmp(inpName,optionNames))
      cfg.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name', inpName)
   end
end
clear pair
clear inpName

% pop up the GUI if every option isn't set -- there are no true defaults from the commandline!
if nArgs/2 ~= length(cfg); 
    
    % NB: The EEGLAB programmers really don't get OOP and so we have to
    % explictly cast our enumeration types, which are subclasses of double,
    % back to double, otherwise we get a weird crash with a cryptic error
    % message
    uilist = { ...
      { 'Style', 'text', 'string', 'Sliding window width (ms)', 'fontweight', 'bold'  } ...
      { 'Style', 'edit', 'string', cfg.window_width_ms, 'tag', 'window_width_ms'} ...
      { 'Style', 'text', 'string', 'Step interval (ms)', 'fontweight', 'bold'  } ...
      { 'Style', 'edit', 'string', cfg.step_width_ms, 'tag', 'step_width_ms'} ...
      ...
      { 'Style', 'text', 'string', 'Decoding dimension', 'fontweight', 'bold'} ...  
      { 'Style', 'popupmenu', 'string', 'spatial|temporal|spatio-temporal' 'tag' 'stmode' 'value' double(cfg.stmode)} ...
      { 'Style', 'text', 'string', 'Averaging mode', 'fontweight', 'bold'  } ...
      { 'Style', 'popupmenu', 'string', 'single-trial|average by condition' 'tag' 'avmode' 'value' double(cfg.avmode)} ...
      ...
      { 'Style', 'text', 'string', 'Analysis mode', 'fontweight', 'bold'  } ...
      { 'Style', 'popupmenu', 'string', 'SVM|SVR' 'tag' 'analysis_mode' 'value' double(cfg.analysis_mode)} ...  
      { 'Style', 'text', 'string', 'Backend', 'fontweight', 'bold'  }  ...
      { 'Style', 'popupmenu', 'string', 'MATLAB|LibSVM|LibLINEAR' 'tag' 'backend' 'value' double(cfg.backend)} ...
      ...
      { 'Style', 'checkbox', 'string' 'Run permutation test' 'tag' 'perm_test' 'value' double(cfg.perm_test)} ...
      { 'Style', 'checkbox', 'string' 'Use RT matching algorithm' 'tag' 'rt_match' 'value' cfg.rt_match} ...
      { 'Style', 'checkbox', 'string' 'Use z-scores' 'tag' 'zscore_convert' 'value' cfg.zscore_convert} ...
      { 'Style', 'checkbox', 'string' 'Extract feature weights' 'tag' 'feat_weights_mode' 'value' cfg.feat_weights_mode} ...
      ... 
      { 'Style', 'text', 'string', 'Cross validation steps', 'fontweight', 'bold'  } ...
      { 'Style', 'edit', 'string', cfg.cross_val_steps, 'tag', 'cross_val_steps'} ...
      { 'Style', 'text', 'string', 'Repetitions thereof', 'fontweight', 'bold'  } ...
      { 'Style', 'edit', 'string', cfg.n_rep_cross_val, 'tag', 'n_rep_cross_val'} ...
      { 'Style', 'text', 'string', 'Permutation repetitions', 'fontweight', 'bold'  } ...
      { 'Style', 'edit', 'string', cfg.permut_rep, 'tag', 'permut_rep'} ...
      ...
      { 'Style', 'checkbox', 'string' 'Show figure' 'tag' 'display_on' 'value' cfg.display_on} ...
      { 'Style', 'checkbox', 'string' 'Include permutation results in figure' 'tag' 'perm_disp' 'value' (cfg.perm_disp & cfg.perm_test) } ...
      { 'Style', 'text', 'string', 'Position zero (ms, for stats and plots)', 'fontweight', 'bold'  } ...
      { 'Style', 'edit', 'string', cfg.pointzero, 'tag', 'pointzero'} ...
    };


    geometry = {...
        [0.4 0.4 0.4 0.4] ...
        [0.4 0.4 0.4 0.4] ...
        [0.4 0.4 0.4 0.4] ...
        [0.2 0.25 0.15 0.2] ... 
        [0.2 0.1 0.2 0.1 0.2 0.1] ... 
        [0.2 0.3 0.2 0.1] ...
    };


    [ outparam userdat strhalt cfg ] = inputgui('geometry', geometry, ...
                                                      'uilist', uilist, ...
                                                      'helpcom', 'pophelp(''pop_dd_decode'');',...
                                                      'title','Run decoding -- pop_dd_decode()');

    if isempty(cfg);
        fprintf('DDTBOX: Nothing to do.\n');
        return;
    end
end

% remember, this is still not exposed publicly
cfg.backend_flags = backend_flags;

cfg.srate = EEG.srate;
cfg.nchannels = EEG.nbchan;
cfg.chaninfo = EEG.chaninfo;
cfg.chanlocs = EEG.chanlocs;

% the GUI functions return strings, so we need to convert just to make sure
% that we have numbers
for fn = fieldnames(cfg)
    key = cell2str(fn);
    value = cfg.(key);
    if ischar(value)
        if strcmp(key,'pointzero') & strcmp(value,'auto')
            % EEG.xmin is in s, so we convert to ms
            % for xmin < 0, this will tell us the offset such that the
            % zeropoint is the time locking event, for xmin > 0, this will 
            % simply give the left edge of the epoch as the zeropoint
            value = max(0, abs(EEG.xmin)) * 1000;
        else 
            value = str2num(value);
        end
        cfg.(key) = value;
    end 
end

% we double check that our numeric values for the different categorical
% options are stored as an enumeration type

cfg.stmode = DDSTMode(cfg.stmode);
cfg.avmode = DDAverage(cfg.avmode);
cfg.analysis_mode = DDMethod(cfg.analysis_mode);
cfg.backend = DDBackend(cfg.backend);


% need to convert numeric fields back to numbers
% need to call decode_erp()
% need to reassemble the pop command
% return the string command
% -------------------------
%com = sprintf('pop_dd_decode( %s, %d, [%s] );', inputname(1), (param3));
%[results, cfg] = decoding_erp(EEG.data, cfg, varargin);
%EEG.ddtbox.results = results;
EEG.ddtbox.cfg = cfg;
% need to fix decode_erp to use new combi structure
% need to fix do_my_classification to use new backend/method distinction
