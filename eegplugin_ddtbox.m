% eegplugin_ddtbox() - EEGLAB plugin for MVPA
%
% Usage:
%   >> eegplugin_ddtbox(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
%
% Notes:
%   This plugins consist of the following Matlab files:
%
% Create a plugin:
%   For more information on how to create an EEGLAB plugin see the
%   help message of eegplugin_besa() or visit http://www.sccn.ucsd.edu/eeglab/contrib.html
%
%
% See also: pop_dd_decode(), pop_dd_analyse(), 
%           pop_dd_display_group(),  pop_dd_display_individual()
%
% Copyright (c) 2016 Phillip Alday and contributors
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


function eegplugin_ddtobox(fig, try_strings, catch_strings)

% create menu
parentmenu  = findobj(fig, 'tag', 'tools');
%eeglabmenu = findobj(fig, 'tag', 'EEGLAB');
%parentmenu = eeglabmenu; 
%helpmenu =  findobj(fig, 'label', 'Help');
%pos = get(helpmenu,'position') - 1;
ddmenu = uimenu(parentmenu, 'label', 'DDTBOX',...
   'separator','on',...
   'tag','DDTBOX', ...
   'userdata', 'startup:on;epoch:on;continuous:on;chanloc:off;');
%set(ddmenu,'position',pos);

uimenu(ddmenu, 'label', 'Single-subject analysis',...
   'separator','off',...,
   'tag','decode', ... 
   'callback',  ...
   [ try_strings.check_epoch_chanlocs '[EEG LASTCOM] = pop_dd_decode(EEG);' catch_strings.add_to_hist ],...
   'userdata', 'epoch:on;continuous:off;chanloc:on;');

uimenu(ddmenu, 'label', 'Single-subject plot',...
   'separator','off',...,
   'tag','display_individual', ... 
   'callback',  ...
   [ try_strings.check_epoch_chanlocs '[EEG LASTCOM] = pop_dd_display_individual(EEG);' catch_strings.add_to_hist ],...
   'userdata', 'epoch:on;continuous:off;chanloc:on;');

uimenu(ddmenu, 'label', 'Group-level analysis',...
   'separator','on',...,
   'tag','analyse', ... 
   'callback',  ...
   [ try_strings.no_check '[EEG LASTCOM] = pop_dd_analyse(EEG);' catch_strings.add_to_hist ],...
   'userdata', 'epoch:on;continuous:off;chanloc:on;');

uimenu(ddmenu, 'label', 'Group-level plot',...
   'separator','off',...,
   'tag','display_group', ... 
   'callback',  ...
   [ try_strings.no_check '[EEG LASTCOM] = pop_dd_display_group(EEG);' catch_strings.add_to_hist ],...
   'userdata', 'epoch:on;continuous:off;chanloc:on;');

