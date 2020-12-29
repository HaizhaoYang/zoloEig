function EigSol_startup()
% EIGSOL_STARTUP  Startup file for EigSol
%   MAKE adds paths of the EigSol to Matlab.

%  Copyright (c) 2016 Yingzhou Li and Haizhao Yang,
%       Stanford University and Duke University 
%  This file is distributed under the terms of the MIT License.

file_path = mfilename('fullpath');
tmp = strfind(file_path,'EigSol_startup');
file_path = file_path(1:(tmp(end)-1));

% Foulder for all soource files recursively
addpath(genpath([file_path 'src']));
addpath(genpath([file_path 'test']));

% Foulder for all external files recursively
addpath(genpath([file_path 'external']));

end