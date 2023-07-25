function FF_out = uigetfile_to_fullpath(varargin)

if nargin>0
[ffname, ffolder] = uigetfile(varargin);
else
    [ffname, ffolder] = uigetfile('*');
end
    
FF_out = fullfile(ffolder,ffname);


end