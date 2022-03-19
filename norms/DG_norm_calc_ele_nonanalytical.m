function [DG_Er]=DG_norm_calc_ele_nonanalytical(u,etpl,ed,etpl_face,coord,E,v)
% Calculating the error in the DG norm
%
% Input(s):
% Dirichlet_Neumann_function - Dirichlet / Neumann function handle
% Dirichlet_function - Dirichlet function handle
% u          - Displacment solution
% etpl       - Element topology structure
% ed         - Degree of freedom steering matrix
% etpl_face  - Element face tolopogy struture (see seed_mesh.m)
% coord      - Coordinates of all nodes in the mesh
% E          - Young's modulus
% v          - Poisson's ratio
%
% Ouput(s):
% DG_Er  - Error for the entire mesh

%  Copyright (C) 2017 Robert Bird
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

% Internal error estimator stiffness computations 
Er1 = DG_strain_area_measure_nonanalytical    (coord,etpl,ed,u,E,v);      
Er3 = DG_calc_jump_in_displacement_flux(coord,etpl,etpl_face,ed,u,E,v);
% -------------------------------------------------------------------------

% Errors associated with the boundary conditions
Er_D   = DG_Dirichlet_BC_nonanalytical(coord,etpl,etpl_face,ed,u,E,v);
Er_DN  = DG_Dirichlet_Neumann_BC_nonanalytical(coord,etpl,etpl_face,ed,u,E,v);
% -------------------------------------------------------------------------

% Summation
Er=Er1+Er_D+Er_DN+Er3;
DG_Er=sqrt(sum(Er));

