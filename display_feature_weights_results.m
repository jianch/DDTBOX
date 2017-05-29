function display_feature_weights_results(ANALYSIS, FW_ANALYSIS)
%
% This function plots results of group-level analyses on feature weights 
% derived using support vector classification or regression.
% 
% This function is called by analyse_feature_weights_erp, but can also be 
% called by custom results plotting scripts.
%
%
% Inputs:
%
%   ANALYSIS         structure containing analysis/plotting settings and
%                    decoding results data
%
%   FW_ANALYSIS      results of the feature weights analyses
%
%
% Outputs:
%
%
% Usage:   display_feature_weights_results(ANALYSIS, FW_ANALYSIS)
%
%
% Copyright (c) 2013-2017 Stefan Bode and contributors
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


%% DISPLAY - MATRIX OF ALL STEPS: Z-STANDARDISED ABSOLUTE FEAUTRE WEIGHTS (FOR DISPLAY)
%__________________________________________________________________________
%
% This matrix is plotted from
% FW_ANALYSIS.AVERAGE_Z_DISP{analysis-time-steps,channel}. Note that this
% can be different from the statistically tested analysis time-windows 

if ANALYSIS.fw.display_matrix == 1

    %create labels
    channel_labels = [];
    for channel = 1:size(FW_ANALYSIS.chanlocs,2)
   
        channel_labels{channel} = FW_ANALYSIS.chanlocs(1,channel).labels;
    
    end
    
    % channels plotted as rows, time-windows as colums
    resorted_data = [];
    resorted_data(:,:) = FW_ANALYSIS.AVERAGE_Z_DISP(:,:);
    resorted_data = resorted_data';

    % create figure
    figure;
    imagesc(resorted_data(:,:));
    hold on;
    
    set(gca,'Ytick',[1:size(FW_ANALYSIS.AVERAGE_Z_DISP,2)]);
    set(gca,'YTickLabel',(channel_labels));
    ylabel('Channel','FontSize',12,'FontWeight','b');
    
    set(gca,'Xtick',[1:size(FW_ANALYSIS.AVERAGE_Z_DISP,1)]);
    set(gca,'XTickLabel',(FW_ANALYSIS.fw_disp));
    xlabel('Analysis time-step','FontSize',12,'FontWeight','b');
    
    title('Z-standardised absolute feature weights','FontSize',14,'FontWeight','b');

end % ANALYSIS.fw.display_matrix


%% DISPLAY - AVERAGE HEAT MAP FOR RELEVANT STEPS: Z-STANDARDISED ABSOLUTE FEAUTRE WEIGHTS (FOR STATS)
%__________________________________________________________________________

if ANALYSIS.fw.display_average_zmap == 1
     
    to_plot = FW_ANALYSIS.AVERAGESTEPS_SELECT_FW_Z_MEAN';
    
    figure;
    topoplot_decoding(to_plot,...
        FW_ANALYSIS.chanlocs,'style','both','electrodes','labelpoint','maplimits','minmax','chaninfo',FW_ANALYSIS.chaninfo, 'colormap', 'jet'); % or: [MIN MAX]
    
    hold on;
    title('Z-standardised absolute feature weights averaged across time-steps','FontSize',10,'FontWeight','b');
    
    clear to_plot;
    
end


%% DISPLAY - AVERAGE HEAT MAP FOR RELEVANT STEPS: THRESHOLD MAP FOR P UNCORRECTED (FOR STATS)
%__________________________________________________________________________

if ANALYSIS.fw.display_average_uncorr_threshmap == 1
    
    to_plot = FW_ANALYSIS.h_matrix_z_averagestep_uncorr;
    
    figure;
    topoplot_decoding(to_plot,...
        FW_ANALYSIS.chanlocs,'style','fill','electrodes','labelpoint','numcontour',1,'conv','off','maplimits',[0 1],'ccolor',[0 0 0],'ecolor',[1 1 1],'chaninfo',FW_ANALYSIS.chaninfo, 'colormap', 'jet');
    
    hold on;
    title('Feature weights uncorrected threshold-map (averaged across time-steps)','FontSize',10,'FontWeight','b');
    
    clear to_plot;
    
end


%% DISPLAY - AVERAGE HEAT MAP FOR RELEVANT STEPS: THRESHOLD MAP FOR P CORRECTED (FOR STATS)
%__________________________________________________________________________

if ANALYSIS.fw.display_average_corr_threshmap == 1
    
    to_plot = FW_ANALYSIS.h_matrix_z_averagestep_corr;
    
    figure;
    topoplot_decoding(to_plot,...
        FW_ANALYSIS.chanlocs,'style','fill','electrodes','labelpoint','numcontour',1,'conv','off','maplimits',[0 1],'ccolor',[0 0 0],'ecolor',[1 1 1],'chaninfo',FW_ANALYSIS.chaninfo, 'colormap', 'jet');
    
    hold on;
    title('Feature weights corrected threshold-map (averaged across time-steps)','FontSize',10,'FontWeight','b');
    
    clear to_plot;
    
end


%% DISPLAY - HEAT MAP FOR EACH RELEVANT STEP: Z-STANDARDISED ABSOLUTE FEATURE WEIGHTS (FOR STATS)
%__________________________________________________________________________

if ANALYSIS.fw.display_all_zmaps == 1
    
    for steps = 1:size(FW_ANALYSIS.p_matrix_z_corr,2)
        
        to_plot = FW_ANALYSIS.AVERAGE_Z_HEATS(steps,:);
        to_plot = to_plot';
        
        figure;
        topoplot_decoding(to_plot,...
            FW_ANALYSIS.chanlocs,'style','both','electrodes','labelpoint','maplimits','minmax','chaninfo',FW_ANALYSIS.chaninfo, 'colormap', 'jet'); % or: [MIN MAX]
        
        hold on;
        title(['Z-standardised absolute feature weights time-step ' num2str(FW_ANALYSIS.fw_analyse(steps))],'FontSize',10,'FontWeight','b');
    
        clear to_plot;
        
    end % step
    
end


%% DISPLAY - HEAT MAP FOR EACH RELEVANT STEP: THRESHOLD MAPS FOR P UNCORRECTED (FOR STATS)
%__________________________________________________________________________

if ANALYSIS.fw.display_all_uncorr_thresh_maps == 1
    
    for steps = 1:size(FW_ANALYSIS.h_matrix_z_uncorr,2)
        
        to_plot(:,:) = FW_ANALYSIS.h_matrix_z_uncorr{steps};
    
        figure;
        topoplot_decoding(to_plot,...
        FW_ANALYSIS.chanlocs,'style','fill','electrodes','labelpoint','numcontour',1,'conv','off','maplimits',[0 1],'ccolor',[0 0 0],'ecolor',[1 1 1],'chaninfo',FW_ANALYSIS.chaninfo, 'colormap', 'jet');
    
        hold on;
        title(['Feature weights uncorrected threshold-map for time-step ' num2str(FW_ANALYSIS.fw_analyse(steps))],'FontSize',10,'FontWeight','b');
    
        clear to_plot;
        
    end % step
    
end


%% DISPLAY - HEAT MAP FOR EACH RELEVANT STEP: THRESHOLD MAPS FOR P CORRECTED (FOR STATS)
%__________________________________________________________________________

if ANALYSIS.fw.display_all_corr_thresh_maps == 1
    
    for steps = 1:size(FW_ANALYSIS.h_matrix_z_corr,2)
        
        to_plot(:,:) = FW_ANALYSIS.h_matrix_z_corr{steps};
    
        figure;
        topoplot_decoding(to_plot,...
        FW_ANALYSIS.chanlocs,'style','fill','electrodes','labelpoint','numcontour',1,'conv','off','maplimits',[0 1],'ccolor',[0 0 0],'ecolor',[1 1 1],'chaninfo',FW_ANALYSIS.chaninfo, 'colormap', 'jet');
    
        hold on;
        title(['Feature weights corrected threshold-map for time-step ' num2str(FW_ANALYSIS.fw_analyse(steps))],'FontSize',10,'FontWeight','b');
    
        clear to_plot;
        
    end % step
    
end

%__________________________________________________________________________