function display_indiv_results_erp(STUDY,RESULTS)
%
% This script gets input from decoding_erp.m and uses specified step-width 
% and results from mutltivariate classification/regression analysis and 
% displays individual results. If permutation test is on and display of 
% permutation results is on, then these results are displayed for comparison.
%
% 
% Inputs:
%
%   STUDY       structure containing participant dataset information and 
%               multivariate classification/regression settings.
%   RESULTS     structure containing decoding results for an individual
%               subject datset.
%
% Outputs:
%
%
% Copyright (c) 2013-2016 Stefan Bode and contributors
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

%% SET GLOBAL VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

global SBJTODO;

%% DISPLAY MAIN RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

% determine the x-axis scaling
% point_zero=floor( ( size(RESULTS.subj_acc,2)/(size(RESULTS.subj_acc,2)*STUDY.step_width * (1000 / SLIST.sampling_rate)) ) * (SLIST.pointzero) );
nsteps = size(RESULTS.subj_acc,2);

for na = 1:size(RESULTS.subj_acc,1)

    figname = ['fig' num2str(na)]; 
    figname = figure('Position',[100 100 800 400]);
    temp_data(1,:) = RESULTS.subj_acc(na,:);
    plot(temp_data,'-ks','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',5);
    hold on;

    xlabel('time-steps [ms]','FontSize',12,'FontWeight','b');
    
    if STUDY.analysis_mode ~= 3
        ylabel('Decoding Accuracy [%]','FontSize',12,'FontWeight','b');
    elseif STUDY.analysis_mode == 3
        ylabel('Fisher-Z correlation coeff','FontSize',12,'FontWeight','b');
    end; 

    XTickLabels(1:1:nsteps) = ( ( (1:1:nsteps) * STUDY.step_width_ms ) - STUDY.step_width_ms) - STUDY.pointzero; 
    point_zero = find(XTickLabels(1,:) == 0);
    
    if STUDY.analysis_mode ~= 3
        line([point_zero point_zero], [100 30],'Color','r','LineWidth',3);
        set(gca,'Ytick',[0:5:100],'Xtick',[1:1:nsteps]);
    elseif STUDY.analysis_mode == 3
        line([point_zero point_zero], [1 -1],'Color','r','LineWidth',3);
        set(gca,'Ytick',[-1:0.2:1],'Xtick',[1:1:nsteps]);
    end
    set(gca,'XTickLabel',XTickLabels);
    
    title(['SBJ' num2str(SBJTODO) ' ' STUDY.dcg_label ' - analysis '...
        num2str(na) ' of ' num2str(size(RESULTS.subj_acc,1))],'FontSize',14,'FontWeight','b');
        
    clear temp_data;

end % na

%% DISPLAY PERMUTATION RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________

if STUDY.perm_disp == 1

    % determine the x-axis scaling
    % point_zero=floor( ( size(RESULTS.subj_perm_acc,2)/(size(RESULTS.subj_perm_acc,2) * STUDY.step_width * (1000 / SLIST.sampling_rate)) ) * (SLIST.pointzero) );
    nsteps = size(RESULTS.subj_perm_acc,2);
    
    for na = 1:size(RESULTS.subj_perm_acc,1)

        figname = ['fig' num2str(na+1)]; 
        figname = figure('Position',[100 100 800 400]);
        temp_data(1,:) = RESULTS.subj_perm_acc(na,:);
        plot(temp_data,'-ks','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',5);
        hold on;

        if STUDY.analysis_mode ~= 3
            ylabel('Decoding Accuracy [%]','FontSize',12,'FontWeight','b');
        elseif STUDY.analysis_mode == 3
            ylabel('Fisher-Z correlation coeff','FontSize',12,'FontWeight','b');
        end; 

        if STUDY.analysis_mode ~= 3
            line([point_zero point_zero], [100 30],'Color','r','LineWidth',3);
            set(gca,'Ytick',[0:5:100],'Xtick',[1:1:nsteps]);
        elseif STUDY.analysis_mode == 3
            line([point_zero point_zero], [1 -1],'Color','r','LineWidth',3);
            set(gca,'Ytick',[-1:0.2:1],'Xtick',[1:1:nsteps]);
        end
        set(gca,'XTickLabel',XTickLabels);

        title(['SBJ' num2str(SBJTODO) ' ' STUDY.dcg_label ' - permutation '...
            num2str(na) ' of ' num2str(size(RESULTS.subj_acc,1))],'FontSize',14,'FontWeight','b');

        clear temp_data;

    end % na

end % perm_disp
%__________________________________________________________________________