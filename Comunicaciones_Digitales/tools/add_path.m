%-----------------------------------------------------------------------------%
%				                FUNDACION FULGOR
%-----------------------------------------------------------------------------%

%-----------------------------------------------------------------------------%
%                                   ADD PATH 
%-----------------------------------------------------------------------------%

proj_dir = mfilename('fullpath');
proj_dir = proj_dir(1: end - length(mfilename));

addpath(genpath([proj_dir, 'TPs/']))
addpath(genpath([proj_dir, 'examples/']))
addpath(genpath([proj_dir, 'sim_tools/']))
addpath(genpath([proj_dir, 'sim_plotter/']))