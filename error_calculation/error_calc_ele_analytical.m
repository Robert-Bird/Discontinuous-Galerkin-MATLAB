function [Er,Er_tot]=error_calc_ele_analytical(u,etpl,ed,etpl_face,coord,E,v)
% Calculating the linear elastic a posteriori error for the mesh
%
% Input(s):
% u          - Displacment solution
% etpl       - Element topology structure
% ed         - Degree of freedom steering matrix
% etpl_face  - Element face tolopogy struture (see seed_mesh.m)
% coord      - Coordinates of all nodes in the mesh
% E          - Young's modulus
% v          - Poisson's ratio
%
% Ouput(s):
% Er     - Error for each active element in the mesh squared
% Er_tot - Error for the entire mesh

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

% Internal error estimator stiffness computations 
Er1 = Error_calc_strong_form_area         (coord,etpl,ed,u,E,v);                
Er2 = Error_calc_jump_in_stress_flux      (coord,etpl,etpl_face,ed,u,E,v);
Er3 = Error_calc_jump_in_displacement_flux(coord,etpl,etpl_face,ed,u,E,v);
% -------------------------------------------------------------------------

% Errors associated with the boundary conditions
Er_D   = Error_calc_Dirichlet_BC(coord,etpl,etpl_face,ed,u,E,v);
Er_DN  = Error_calc_Dirichlet_Neumann_BC(coord,etpl,etpl_face,ed,u,E,v);
Er_N   = Error_calc_Neumann_BC_analytical(coord,etpl,etpl_face,ed,u,E,v);
% -------------------------------------------------------------------------

% Summation
Er=Er1+Er2+Er3+Er_D+Er_DN+Er_N;
Er_tot=sqrt(sum(Er));

