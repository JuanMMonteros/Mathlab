%-------------------------------------------------------------------------%
%				              FUNDACION FULGOR
%-------------------------------------------------------------------------%
%                                  ADD PATH 
%-------------------------------------------------------------------------%

proj_dir = mfilename('fullpath');
proj_dir = proj_dir(1: end - length(mfilename));

addpath(genpath([proj_dir, 'Module_test/']))
addpath(genpath([proj_dir, 'Simu_test/']))
addpath(genpath([proj_dir, 'SRC/']))
addpath(genpath([proj_dir, 'tests/']))
addpath(genpath([proj_dir, 'Tools/']))