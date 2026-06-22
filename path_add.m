% path_add: Adds all subfolders to the solution path
% Author: Robert Bird
% Date:   11/06/2017
% Description: Adds paths and all subpaths to MATLAB:
% adaptivity/
% boundary_conditions/
% error_calculation/
% mesh_generation/
% mesh_generation/square/
% Mesh2d/ - Comment in line 38 if Mesh2d is installed
% norms/
% plotters/
% rhs/
% solvers/
% stiffness_matrix/
% stiffness_matrix/shape_func/
% stiffness_matrix/shape_func/shape_function_generation/
% stiffness_matrix/volume_matrix/
% stiffness_matrix/surface_matrix/
% problem_set_up/
% tests/
% example_problems/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function path_add
addpath('adaptivity');
addpath('boundary_conditions'); 
addpath('error_calculation'); 
addpath('mesh_generation');
addpath('mesh_generation/square');
%addpath(genpath('Mesh2d'));
addpath('norms/');
addpath('plotters');
addpath('rhs');
addpath('solvers'); 
addpath('stiffness_matrix'); 
addpath('stiffness_matrix/shape_func/');
addpath('stiffness_matrix/shape_func/shape_function_generation');
addpath('stiffness_matrix/volume_matrix');
addpath('stiffness_matrix/surface_matrix');
addpath(genpath('problem_set_up'));
addpath(genpath('tests'));
addpath(genpath('example_problems'));
