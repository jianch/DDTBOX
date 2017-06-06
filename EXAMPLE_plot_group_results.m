% EXAMPLE_plot_group_results.m
%
% This script is used for plotting group-level results of 
% classifier accuracy/SVR performance and feature weights analyses.
% Group-level statistical analyses must be performed before running this script.
% Please make copies of this script for your own projects.
%
% This script calls display_group_results_erp and
% display_feature_weights_results, both of which are included in DDTBOX.
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


%% General settings

% Full filepath of group-level results file
group_results_file = ['/Users/danielfeuerriegel/Desktop/DDTBOX Project/MVPA_WORKSHOP/DECODING_RESULTS/Tutorial/FLANKER_GROUPRES_NSBJ5_win10_steps10_av1_st3_SVM_LIBSVM_DCGCorrect vs. Error.mat'];

% Load the data file to get ANALYSIS parameters
load(group_results_file);

% Select whether to plot decoding performance and feature weights results
PLOT.decoding_performance_results = 1; % 1 = plot; 0 = don't plot
PLOT.feature_weights_results = 1; % 1 = plot; 0 = don't plot


%% Settings for group decoding performance results (classifier accuracy/SVR performance)

% From the group-level analyses configuration script:
ANALYSIS.permdisp = 1; % display the results from permutation test in figure as separate line? 0=no / 1=yes
ANALYSIS.disp.sign = 1; % display statistically significant steps in results figure? 0=no / 1=yes
ANALYSIS.plot_robust = 0; % Choose estimate of location to plot. 0 = arithmetic mean / 1 = trimmed mean / 2 = median
% Note: You can plot the mean even if the data were originally plotted with
% the median or trimmed mean. If you originally plotted the trimmed mean,
% then you can also plot the median in this script.
ANALYSIS.disp.temporal_decoding_colormap = 'jet'; % Colormap for temporal decoding results scalp maps

% figure position on the screen
PLOT.FigPos = [100 100 800 400];

% Figure title settings
PLOT.TitleFontSize = 14;
PLOT.TitleFontWeight = 'Bold'; % 'Normal' (Regular) or 'Bold'

% Title text automatically generated based on the decoding method used
if ANALYSIS.stmode == 1 && ANALYSIS.analysis_mode ~=3 % Spatial SVM Classification
    PLOT.TitleString = 'Spatial SVM ';
elseif ANALYSIS.stmode == 2 && ANALYSIS.analysis_mode ~=3 % Temporal SVM Classification
    PLOT.TitleString = 'Temporal SVM ';
elseif ANALYSIS.stmode == 3 && ANALYSIS.analysis_mode ~=3 % Spatiotemporal SVM Classification
    PLOT.TitleString = 'Spatiotemporal SVM ';
elseif ANALYSIS.stmode == 1 && ANALYSIS.analysis_mode ==3 % Spatial SVR
    PLOT.TitleString = 'Spatial SVR';  
elseif ANALYSIS.stmode == 2 && ANALYSIS.analysis_mode ==3 % Temporal SVR
    PLOT.TitleString = 'Temporal SVR '; 
elseif ANALYSIS.stmode == 3 && ANALYSIS.analysis_mode ==3 % Spatiotemporal SVR
    PLOT.TitleString = 'Spatiotemporal SVR ';  
end % of if ANALYSIS.stmode


%% Decoding Performance X and Y Axis Properties

% Axis label properties
PLOT.xlabel.FontSize = 12;
PLOT.ylabel.FontSize = 12;
PLOT.xlabel.FontWeight = 'Bold'; % 'Normal' (Regular) or 'b' / 'Bold'
PLOT.ylabel.FontWeight = 'Bold'; % 'Normal' (Regular) or 'b' / 'Bold'

% X axis label text
PLOT.xlabel.Text = 'Time-steps [ms]';

% Y axis label text
if ANALYSIS.analysis_mode ~= 3 % If not using SVR
    
    PLOT.ylabel.Text = 'Classification Accuracy [%]';
    
elseif ANALYSIS.analysis_mode == 3 % If using SVR
    
    PLOT.ylabel.Text = 'Fisher-transformed correlation coeff';
    
end % of if ANALYSIS.analysis_mode

% _________________________________________________________________________

% Define x/y-axis limits and tick marks

% Y-axis depends on analysis mode
if ANALYSIS.analysis_mode ~= 3 % If not using SVR
    
    PLOT.Y_min = 40; % Y axis lower bound (in % accuracy)
    PLOT.Y_max = 80; % Y axis upper bound (in % accuracy)
    PLOT.Ysteps = 5; % Interval between Y axis labels/tick marks
    
elseif ANALYSIS.analysis_mode == 3 % If using SVR
    
    PLOT.Y_min = -0.5; % Y axis lower bound (Fisher-Z corr coeff)
    PLOT.Y_max = 0.5; % Y axis upper bound (Fisher-Z corr coeff)
    PLOT.Ysteps = 0.1; % Interval between Y axis labels/tick marks
    
end % of if ANALYSIS.analysis_mode

PLOT.X_min = 1; % X axis lower bound (first time point)
PLOT.X_max = ANALYSIS.xaxis_scale(2,end); % Maximum value of X axis value. 
% Default is ANALYSIS.xaxis_scale(2,end) which is the last time window selected for group-level analyses.

% Automated calculations (no input required)
PLOT.Xsteps = ANALYSIS.step_width_ms;
PLOT.Ytick = [PLOT.Y_min:PLOT.Ysteps:PLOT.Y_max];
PLOT.Xtick = [ANALYSIS.xaxis_scale(1,1) : ANALYSIS.xaxis_scale(1,end)];
PLOT.XtickLabel = ANALYSIS.xaxis_scale(2,:) - ANALYSIS.pointzero; 


%% Define properties of lines showing decoding performance and error bars

% Actual decoding results
PLOT.Res.Line = '-ks'; % Line colour and style
PLOT.Res.LineWidth = 2;
PLOT.Res.MarkerEdgeColor = 'k';
PLOT.Res.MarkerFaceColor = 'w';
PLOT.Res.MarkerSize = 5;

% Error bar plotting
PLOT.Res.Error = 'k'; % Line colour and style
PLOT.Res.ErrorLineWidth = 0.5;
PLOT.Res.ErrorLine = 'none'; % Disables lines between error bars across steps

% Properties of line showing permutation / chance results
PLOT.PermRes.Line = '-ks'; % Line colour and style
PLOT.PermRes.LineWidth = 2;
PLOT.PermRes.MarkerEdgeColor = 'b';
PLOT.PermRes.MarkerFaceColor = 'w';
PLOT.PermRes.MarkerSize = 5;

% Error bars for chance/permutation results
PLOT.PermRes.Error = 'b'; % Line colour and style
PLOT.PermRes.ErrorLineWidth = 0.5;
PLOT.PermRes.ErrorLine = 'none'; % Disables lines between error bars across steps


%% Decoding performance plot annotations

% Define properties of line showing event onset
PLOT.PointZero.Color = 'r'; % Colour of line denoting event onset
PLOT.PointZero.LineWidth = 3; % Width of line denoting event onset
PLOT.PointZero.Point = find(ANALYSIS.data(3,:) == 1);

% Define properties of statistical significance markers
PLOT.Sign.LineColor = 'y';
PLOT.Sign.LineWidth = 10;

% Positions of statistical significance markers
if ANALYSIS.analysis_mode ~= 3 % If not using SVR
    PLOT.Sign.LinePos = [PLOT.Y_min + 0.5, PLOT.Y_max - 0.5];
elseif ANALYSIS.analysis_mode == 3 % If using SVR
    PLOT.Sign.LinePos = [PLOT.Y_min, PLOT.Y_max];
end % of if ANALYSIS.analysis_mode


%% Plotting settings for feature weights analyses

% Maps and stats for averaged analysis time windows
ANALYSIS.fw.display_average_zmap = 1; % z-standardised average FWs
ANALYSIS.fw.display_average_uncorr_threshmap = 1; % thresholded map uncorrected t-test results
ANALYSIS.fw.display_average_corr_threshmap = 1; % thresholded map t-test results corrected for multiple comparisons

% Maps and stats for each analysis time window
ANALYSIS.fw.display_all_zmaps = 1; % z-standardised average FWs
ANALYSIS.fw.display_all_uncorr_thresh_maps = 0; % thresholded map uncorrected t-test results
ANALYSIS.fw.display_all_corr_thresh_maps = 0; % thresholded map t-test results corrected for multiple comparisons

% Extra plotting options:
ANALYSIS.fw.colormap = 'jet'; % Colormap for plotting of scalp feature weights maps



%% Plot the classification accuracy/SVR performance results
if PLOT.decoding_performance_results
    display_group_results_erp(ANALYSIS, PLOT);
end % of if PLOT.decoding_performance

%% Plot feature weights results
if PLOT.feature_weights_results
     display_feature_weights_results(ANALYSIS, FW_ANALYSIS);
end % of if PLOT.feature_weights


