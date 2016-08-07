% DDDrawmode.m
%
% This file provides an enumeration for the different resampling modes
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

classdef DDDrawmode < DDEnum
    enumeration
        % Testing against: 
        % average permutated distribution (faster)
        average (1)
        % random values drawn form permuted distribution (stricter)
        random (2)
    end  
end