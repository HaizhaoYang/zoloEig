function MuFiM_startup()
% MUFIM_STARTUP  Startup file for MuFiM
%   MUFIM_STARTUP() adds paths of the MultiFrontal to Matlab.

%   Copyright (c) 2016 Yingzhou Li, Stanford University

file_path = mfilename('fullpath');
tmp = strfind(file_path,'MuFiM_startup');
file_path = file_path(1:(tmp(end)-1));

addpath(genpath([file_path 'src']));

addpath(genpath([file_path 'external']));

end
