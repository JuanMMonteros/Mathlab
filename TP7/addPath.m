%-------------------------------------------------------------------------%
%				              FUNDACION FULGOR
%-------------------------------------------------------------------------%
%                                  ADD PATH 
%-------------------------------------------------------------------------%

proj_dir = mfilename('fullpath');
proj_dir = proj_dir(1: end - length(mfilename));

addpath(genpath([proj_dir, 'cpr_and_bps/']))
addpath(genpath([proj_dir, 'EJ1/']))
addpath(genpath([proj_dir, 'EJ2/']))
addpath(genpath([proj_dir, 'EJ3/']))