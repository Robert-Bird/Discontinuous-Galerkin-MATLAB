%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Robert Bird
% Date:   11/06/2017
% Description: Adds paths:
% -/Adaptivity
% -/Boundary conditions
% -/Deferinement
% -/Error Calculation
% -/Eshelby stress
% -/LHS
% -/Mesh_Refine
% -/norms
% -/Plotters
% -/Propagation
% -/Solvers
% -/Stiffness matrix
% -/Stiffness matrix/Shape_func
% -/triangle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function path_add
addpath('adaptivity');
addpath('boundary_conditions'); 
addpath('error_calculation'); 
addpath('mesh_generation');
addpath('mesh_generation/square');
addpath(genpath('Mesh2d'));
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
