% EXAMPLE_analyse_decoding_results.m
%
% This script is used for specifying configuration settings for analyses of
% decoding results using DDTBOX. All analysis parameters are specified here
% and passed to analyse_decoding_erp.m
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

%% Housekeeping
clear variables;
close all;


%% Filepaths and locations of subject datasets

% Enter the name of the study (for labeling saved decoding results files)
study_name = 'EXAMPLE';

% Base directory path (where single subject EEG datasets and channel locations information are stored)
bdir = '/Users/danielfeuerriegel/Desktop/DDTBOX Project/MVPA_WORKSHOP/';

% Output directory (where decoding results have been saved)
output_dir = '/Users/danielfeuerriegel/Desktop/DDTBOX Project/MVPA_WORKSHOP/DECODING_RESULTS/Test/';

%% Select subject datasets and discrimination groups

% Set which subjects datasets to decode
sbjs_todo = [1:4];

% Enter the discrimination group for classification. Two discrimination
% groups can be entered when using cross-condition decoding.
dcg_todo = [1];

% Perform cross-condition decoding? 1 = yes / 0 = no
cross = 0;

%% EEG dataset information

nchannels = 64; % number of channels
channel_names_file = 'channel_inf.mat'; % Name of .mat file containing channel labels and channel locations
channellocs = [bdir, 'locations/']; % Path of directory containing channel information file
sampling_rate = 1000; % Data sampling rate in Hz
pointzero = 100; % Corresponds to the time of the event/trigger code relative to the prestimulus baseline (in ms)


%% Condition and discrimination group (dcg) information

% Label each condition
% Usage: cond_labels{condition number} = 'Name of condition';
% Example: cond_labels{1} = 'Correct Responses';
% Condition label {X} corresponds to data in column X of the single subject
% data arrays.
cond_labels{1} = 'condition_A';
cond_labels{2} = 'condition_B';
cond_labels{3} = 'condition_C';
cond_labels{4} = 'condition_D';
        
% Discrimination groups
% Enter the condition numbers of the conditions to discriminate between
% Usage: dcg{discrimination group number} = [condition 1, condition 2];
% Example: dcg{1} = [1, 2]; to compare conditions 1 and 2 for dcg 1
dcg{1} = [1, 3]; 
dcg{2} = [2, 4]; 
              
% Label each discrimination group
% Usage: dcg_labels{Discrimination group number} = 'Name of discrimination group'
% Example: dcg_labels{1} = 'Correct vs. Error Responses';
dcg_labels{1} = 'A vs. C';
dcg_labels{2} = 'B vs. D';


% This section automaticallly fills in various parameters related to dcgs and conditions 
ndcg = size(dcg, 2);
nclasses = size(dcg{1}, 2);      
ncond = size(cond_labels, 2);



%% Group-level analysis and plotting parameters

% define all parameters of results to analyse & Plot
%______________________________________________________________________
allchan = 1; % Are all possible channels analysed? 1=yes (default if spatial/spatio-temporal) / 2=no
relchan = []; % specify channels to be analysed (for temporal only)

analysis_mode = 1; % ANALYSIS mode (1=SVM with LIBSVM / 2=SVM with liblinear / 3=SVR with LIBSVM)
stmode = 1; % SPACETIME mode (1=spatial / 2=temporal / 3=spatio-temporal)
avmode = 1; % AVERAGE mode (1=no averaging; single-trial / 2=run average) 
window_width_ms = 50; % width of sliding window in ms
step_width_ms = 50; % step size with which sliding window is moved through the trial

pstats = 0.05; % critical p-value
group_level_analysis = 2; % Select statistical analysis method: 1 = Global null and prevalence testing based on the minimum statistic / 2 = Global null testing with t tests

% If using minimum statistic approach for group-level analyses:
P2 = 100000; % Number of second-level permutations to use
minstat_multcomp = 1; % Correct for multiple comparisons using the maximum statistic approach:
% 0 = no correction
% 1 = correction based on the maximum statistic (also applied to prevalence lower bound testing)

% If using t test approach for group-level analyses:
permstats = 2; % Testing against: 1=theoretical chance level / 2=permutation test results
drawmode = 1; % Testing against: 1=average permutated distribution (default) / 2=random values drawn form permuted distribution (stricter)
groupstats_ttest_tail = 'right'; % Choose between two-tailed or one-tailed tests. 'both' = two-tailed / 'right' / 'left' = one-tailed testing for above/below chance accuracy
use_robust = 0; % Use Yuen's t, the robust version of the t test? 1 = Yes / 0 = No
trimming = 20; % If using Yuen's t, select the trimming percentage for the trimmed mean (20% recommended)

multcompstats = 0; % Correction for multiple comparisons: 
% 0 = no correction
% 1 = Bonferroni correction
% 2 = Holm-Bonferroni correction
% 3 = Strong FWER Control Permutation Test
% 4 = Cluster-Based Permutation Test
% 5 = KTMS Generalised FWER Control
% 6 = Benjamini-Hochberg FDR Control
% 7 = Benjamini-Krieger-Yekutieli FDR Control
% 8 = Benjamini-Yekutieli FDR Control
n_iterations = 1000; % Number of permutation or bootstrap iterations for resampling-based multiple comparisons correction procedures
ktms_u = 2; % u parameter of the KTMS GFWER control procedure
cluster_test_alpha = 0.05; % For cluster-based test: Significance threshold for detecting effects at individual time windows (e.g. 0.05)

% Group-level classifier accuracy results plotting options
disp.on = 1; % display a results figure? 0=no / 1=yes
permdisp = 1; % display the results from permutation test in figure as separate line? 0=no / 1=yes
disp.sign = 1; % display statistically significant steps in results figure? 0=no / 1=yes
plot_robust = 0; % Choose estimate of location to plot. 0 = arithmetic mean / 1 = trimmed mean / 2 = median
plot_robust_trimming = 20; % Percent to trim if using the trimmed mean

% Feature weight analysis options
fw.do = 0; % analyse feature weights? 0=no / 1=yes
fw.corrected = 1; % Use feature weights corrected using Haufe et al. (2014) method? 0=no / 1=yes
use_robust_fw = 0; % Use Yuen's t, the robust version of the t test for feature weights? 1 = Yes / 0 = No
trimming_fw = 20; % If using Yuen's t, select the trimming percentage for the trimmed mean (20% recommended)
fw_ttest_tail = 'right';
fw.multcompstats = 1; % Feature weights correction for multiple comparisons:
% 1 = Bonferroni correction
% 2 = Holm-Bonferroni correction
% 3 = Strong FWER Control Permutation Test
% 4 = Cluster-Based Permutation Test (Currently not available)
% 5 = KTMS Generalised FWER Control
% 6 = Benjamini-Hochberg FDR Control
% 7 = Benjamini-Krieger-Yekutieli FDR Control
% 8 = Benjamini-Yekutieli FDR Control

% if feature weights are analysed, specify what is displayed
%__________________________________________________________________

% 0=no / 1=yes
fw.display_matrix = 0; % feature weights matrix

% maps and stats for averaged analysis time windows
fw.display_average_zmap = 0; % z-standardised average FWs
fw.display_average_uncorr_threshmap = 0; % thresholded map uncorrected t-test results
fw.display_average_corr_threshmap = 0; % thresholded map t-test results corrected for multiple comparisons

% maps and stats for each analysis time window
fw.display_all_zmaps = 0; % z-standardised average FWs
fw.display_all_uncorr_thresh_maps = 0; % thresholded map uncorrected t-test results
fw.display_all_corr_thresh_maps = 0; % thresholded map t-test results corrected for multiple comparisons



%% Copy all settings into a structure
% This structure is passed as a single input argument to
% analyse_decoding_erp

% DF TODO: Move this section into a separate script, however kept in this
% script for now to assist debugging.

ANALYSIS.bdir = bdir;
ANALYSIS.output_dir = output_dir;
ANALYSIS.nchannels = nchannels;
ANALYSIS.channel_names_file = channel_names_file;
ANALYSIS.channellocs = channellocs;
ANALYSIS.sampling_rate = sampling_rate;
ANALYSIS.pointzero = pointzero;
ANALYSIS.cond_labels = cond_labels;
ANALYSIS.dcg = dcg;
ANALYSIS.dcg_labels = dcg_labels;
ANALYSIS.ndcg = ndcg;
ANALYSIS.nclasses = nclasses;
ANALYSIS.ncond = ncond;
ANALYSIS.study_name = study_name;
ANALYSIS.sbjs_todo = sbjs_todo;
ANALYSIS.dcg_todo = dcg_todo;
ANALYSIS.cross = cross;
ANALYSIS.allchan = allchan;
ANALYSIS.relchan = relchan;
ANALYSIS.analysis_mode = analysis_mode;
ANALYSIS.stmode = stmode;
ANALYSIS.avmode = avmode;
ANALYSIS.window_width_ms = window_width_ms;
ANALYSIS.step_width_ms = step_width_ms;
ANALYSIS.pstats = pstats;
ANALYSIS.group_level_analysis = group_level_analysis;
ANALYSIS.P2 = P2;
ANALYSIS.minstat_multcomp = minstat_multcomp;
ANALYSIS.permstats = permstats;
ANALYSIS.drawmode = drawmode;
ANALYSIS.groupstats_ttest_tail = groupstats_ttest_tail;
ANALYSIS.use_robust = use_robust;
ANALYSIS.trimming = trimming;
ANALYSIS.multcompstats = multcompstats;
ANALYSIS.n_iterations = n_iterations;
ANALYSIS.ktms_u = ktms_u;
ANALYSIS.cluster_test_alpha = cluster_test_alpha;
ANALYSIS.disp.on = disp.on;
ANALYSIS.permdisp = permdisp;
ANALYSIS.disp.sign = disp.sign;
ANALYSIS.plot_robust = plot_robust;
ANALYSIS.plot_robust_trimming = plot_robust_trimming;
ANALYSIS.fw.do = fw.do;
ANALYSIS.fw.corrected = fw.corrected;
ANALYSIS.use_robust_fw = use_robust_fw;
ANALYSIS.trimming_fw = trimming_fw;
ANALYSIS.fw_ttest_tail = fw_ttest_tail;
ANALYSIS.fw.multcompstats = fw.multcompstats;
ANALYSIS.fw.display_matrix = fw.display_matrix;
ANALYSIS.fw.display_average_zmap = fw.display_average_zmap;
ANALYSIS.fw.display_average_uncorr_threshmap = fw.display_average_uncorr_threshmap;
ANALYSIS.fw.display_average_corr_threshmap = fw.display_average_corr_threshmap;
ANALYSIS.fw.display_all_zmaps = fw.display_all_zmaps;
ANALYSIS.fw.display_all_uncorr_thresh_maps = fw.display_all_uncorr_thresh_maps;
ANALYSIS.fw.display_all_corr_thresh_maps = fw.display_all_corr_thresh_maps;


%% Analyse decoding results for specified subjects and dcgs

analyse_decoding_erp(ANALYSIS);

