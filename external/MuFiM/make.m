function make(installflag)
% MAKE  Make file for MuFiM
%   MAKE, MAKE(1) install MuFiM under the userpath/MuFiM. The path to
%   installed MuFiM will be added in 'startup.m' such that the MuFiM will
%   be included automatically.
%
%   MAKE(-1) uninstall MuFiM.
%

%   Copyright (c) 2016 Yingzhou Li, Stanford University

if( nargin == 0 )
    installflag = 1;
end

src_path = ['src' filesep];
install_path = strsplit(userpath,pathsep);
matlab_startup_file = [install_path{1} filesep 'startup.m'];
install_path = [install_path{1} filesep 'MuFiM' filesep];
install_src_path = [install_path filesep 'src' filesep];

% make clear, installflag<0
if(installflag < 0)
    if(exist(install_path, 'dir'))
        rmpath([install_path 'src' filesep]);
        rmdir(install_path,'s');
    end
    if(exist(matlab_startup_file, 'file'))
        startup_path_data = importdata(matlab_startup_file);
        fid = fopen(matlab_startup_file, 'w+');
        for i=1:length(startup_path_data)
            if(~strcmp(startup_path_data{i},['run ' install_path 'MuFiM_startup.m']))
                fprintf(fid,'%s\n',startup_path_data{i});
            end
        end
        fclose(fid);
    end
    return;
end

% make install

if(installflag == 1)
    if(~exist(install_path, 'dir'))
        mkdir(install_path);
    end
    if(~exist(install_src_path, 'dir'))
        mkdir(install_src_path);
    end
    copyfile('MuFiM_startup.m',install_path);
    copyfile([src_path '*'],install_src_path);
    
    flagexist = 0;
    if exist(matlab_startup_file,'file')
        startup_path_data = importdata(matlab_startup_file);
        for i=1:length(startup_path_data)
            if(strcmp(startup_path_data{i},...
                ['run ' install_path 'MuFiM_startup.m']))
                flagexist = 1;
            end
        end
    end
        
    if(~flagexist)
        fid = fopen(matlab_startup_file, 'at');
        fprintf(fid, '%s\n', ['run ' install_path 'MuFiM_startup.m']);
        fclose(fid);
    end
    run([install_path 'MuFiM_startup.m']);
    return;
end

end
