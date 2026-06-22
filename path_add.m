% path_add: Adds all subfolders to the solution path
% Author: Robert Bird
% Date:   11/06/2017
% Description: Adds paths and all subpaths to MATLAB:
% adaptivity/
% boundary_conditions/
% error_calculation/
% mesh_generation/
% mesh_generation/square/
% Mesh2d/ - third-party (separate licence); auto-detected if added (see end)
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

% Mesh2D is third-party (separate licence) and is NOT shipped with this
% project. If you have added a Mesh2d folder into the main directory, put it
% (and its subfolders) on the path and install the size-compatibility shim
% into Mesh2d/private so refine2 works on recent MATLAB. This affects Mesh2D
% only; nothing here modifies a Mesh2D source file.
if isfolder('Mesh2d')
    ws = warning('off','MATLAB:dispatcher:nameConflict');   % quiet the install-time shadow notice only
    addpath(genpath('Mesh2d'));
    install_mesh2d_shim('Mesh2d');
    warning(ws);                                            % restore previous warning state
end
end
