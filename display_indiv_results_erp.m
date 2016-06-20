function display_indiv_results_erp(STUDY,RESULTS)
%__________________________________________________________________________
% DDTBOX script written by Stefan Bode 01/03/2013
%
% The toolbox was written with contributions from:
% Daniel Bennett, Daniel Feuerriegel, Phillip Alday
%
% The author further acknowledges helpful conceptual input/work from: 
% Jutta Stahl, Simon Lilburn, Philip L. Smith, Elaine Corbett, Carsten Murawski, 
% Carsten Bogler, John-Dylan Haynes
%__________________________________________________________________________
%
% Gets input from DECODING_ERP.m
% Uses specified step-width and results from mutltivariate classification/
% regression analysis and displays individual results.
% If premutation test is on and display of permutation results is on, then
% these results are displayed for comparison.

%__________________________________________________________________________
%
% Variable naming convention: STRUCTURE_NAME.example_variable

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