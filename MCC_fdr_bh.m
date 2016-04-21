function [fdr_corrected_h] = MCC_fdr_bh(p_values, alpha_level)

%__________________________________________________________________________
% Multiple comparisons correction function written by Daniel Feuerriegel 21/04/2016 
% to complement DDTBOX scripts written by Stefan Bode 01/03/2013.
%
% The toolbox was written with contributions from:
% Daniel Bennett, Jutta Stahl, Daniel Feuerriegel, Phillip Alday
%
% The author (Stefan Bode) further acknowledges helpful conceptual input/work from: 
% Simon Lilburn, Philip L. Smith, Elaine Corbett, Carsten Murawski, 
% Carsten Bogler, John-Dylan Haynes
%__________________________________________________________________________
%
% This script receives a vector of p-values and outputs
% false discovery rate corrected null hypothesis test results (Benjamin-Hochberg procedure).
% The number of tests is determined by the length of the vector of p-values.
%
%
% requires:
% - p_values (vector of p-values from the hypothesis tests of interest)
% - alpha_level (uncorrected alpha level for statistical significance)
%
%
% outputs:
% fdr_corrected_h (vector of false discovery rate corrected hypothesis tests 
% derived from comparing p-values to false discovery rate adjusted critical alpha level. 
% 1 = statistically significant, 0 = not statistically significant)
%__________________________________________________________________________
%
% Variable naming convention: STRUCTURE_NAME.example_variable

n_total_comparisons = length(p_values); % Get the number of comparisons
fdr_corrected_h = zeros(1, length(p_values)); % preallocate

sorted_p = sort(p_values); % Sort p-values from smallest to largest

% Find critical k value
for benhoch_step = 1:n_total_comparisons
    if sorted_p(benhoch_step) <= (benhoch_step / n_total_comparisons) * alpha_level
        benhoch_critical_alpha = sorted_p(benhoch_step);
    end
end

% If no steps are significant set critical alpha to zero
if ~exist('benhoch_critical_alpha', 'var')
    benhoch_critical_alpha = 0;
end

% Declare tests significant if they are smaller than or equal to the adjusted critical alpha
fdr_corrected_h(p_values <= benhoch_critical_alpha) = 1;