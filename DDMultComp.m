% DDMultComp.m
%
% This file provides an enumeration for the different multiple-comparisons
% corrections
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

classdef DDMultComp < DDEnum
    enumeration
        % no correction
        none (1)
        % Bonferroni correction
        bonferroni (2)
        % Holm-Bonferroni correction
        holm_bonferroni (3)
        % Strong FWER Control Permutation Test
        fwer (4) 
        % Cluster-Based Permutation Test
        cluster (5) 
        % KTMS Generalised FWER Control
        ktms (6) 
        % Benjamini-Hochberg FDR Control
        hochberg_fdr (7) 
        % Benjamini-Krieger-Yekutieli FDR Control
        krieger_fdr (8)
        % Benjamini-Yekutieli FDR Control
        yekutieli_fdr (9)
    end
end

