function [K,tndof,ed]=DG_algorithm(etpl,etpl_face,coord,E,v)
% Stiffness matrix generation
%
% Input(s):
% etpl       - Element tolopogy struture 
% etpl_face  - Element face tolopogy struture 
% coord      - Element coordinates 
% E          - Youngs modulus
% v          - Poisson's ratio
%
% Ouput(s):
% K                - Sparse global stiffness matrix
% tndof            - Total number of degrees of freedom
% ed               - Degree of freedom steering matrix

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

[kval,krow,kcol,ed] = vol_int(coord,etpl,E,v);                             % Local volumetric stiffness calculation
[k] = surf_int_new(coord,etpl,etpl_face,E,v,ed);                           % Local internal stiffness calculation
[k_D]=surf_dirichlet(coord,etpl,etpl_face,E,v,ed);                         % Local Dirichlet BC external stiffness calculation
[k_DN]=surf_dirichlet_newmann(coord,etpl,etpl_face,E,v,ed);                % Local Dirichlet Newmann BC external stiffness calculation
K = sparse([krow;k.r1;k.r2;k.r3;k.r4;k_D.r1;k_DN.r1],...                   % Global sparse stiffness matrix storage
    [kcol;k.c1;k.c2;k.c3;k.c4;k_D.c1;k_DN.c1],...                          %
    [kval;k.v1;k.v2;k.v3;k.v4;k_D.v1;k_DN.v1],max(ed(:)),max(ed(:)));      %

tndof=max(ed(:));                                                          % total number of degrees of freedom
