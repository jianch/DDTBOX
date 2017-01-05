% DDMultComp.m
%
% This file provides a covenient base class for enumerations
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

classdef DDEnum < double

    methods
        function e = DDEnum(i)
            if isa(i,'char')
               i =  str2num(i);
            end
            e@double(i)
        end
        
        function s = name(obj)
            s = [class(obj) '.' char(obj)];
        end
        
        function s = val(obj)
            s = char(obj);
        end
        
        function disp(obj)
            disp( name(obj) )
        end
    end
    
    methods(Static)
        function s = enum2str(obj,delimiter)
            [~, names] = enumeration(obj);
            % join() was introduced in 2013b
            s = join(names,delimiter);
        end

    end
    
end

