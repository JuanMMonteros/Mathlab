%-----------------------------------------------------------------------------%
%				                FUNDACION FULGOR
%-----------------------------------------------------------------------------%

%-----------------------------------------------------------------------------%
%                                   ADD PATH 
%-----------------------------------------------------------------------------%

proj_dir = mfilename('fullpath');
proj_dir = proj_dir(1: end - length(mfilename));

addpath(genpath([proj_dir, 'module_tests/']))
addpath(genpath([proj_dir, 'sim_tests/']))
addpath(genpath([proj_dir, 'src/']))
addpath(genpath([proj_dir, 'tests/']))
addpath(genpath([proj_dir, 'tools/']))