% This provides a compatibility wrapper to address the shadowing of the 
% Matlab included strjoin() function by BCILAB and CleanLine (which
% includes parts of BCILAB) 

function result = ddwrap_strjoin(C, delimiter)
% Try the builtin order for arguments, if that fails try the BCILAB order
% Note that from MATLAB 2016b, you can simply use join().
% Although join() exists from 2013b, it wasn't overloaded to replace
% strjoin() until 2016b.

    try
        result = strjoin(C, delimiter);
    catch 
        result = strjoin(delimiter, C);
    end
end