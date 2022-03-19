function [coord,etpl,etpl_face]=hp_adapt(Er_ele,etpl,etpl_face,coord,d_max,d_min)
% hp-adaptivity applied to the mesh.
%
% Input(s):
% Er_ele    - Local error for all elements in the mesh
% etpl      - Element topology structure
% etpl_face - Element face topology matrix
% coord     - Element coordinates
% d_max     - d_1 for hp adaptive scheme
% d_min     - d_2 for hp adaptive scheme
%
% Ouput(s):
% coord     - New element coordinates
% etpl      - New element topology structure
% etpl_face - New element face topology matrix
%
%  See also H_ADAPT, P_ADAPT.

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

% Creation of flags for p and h adaptivity     
[refine_flag]=h_adapt(Er_ele,etpl,etpl_face,d_max);                        % h-flag adaptive function
etpl=p_adapt(Er_ele,etpl,d_max,d_min);                                     % p-flag adaptive function

%refine_flag(:,2)=etpl.tree(:,1);
[coord,etpl,etpl_face]=mesh_refine_topology_marking(coord,etpl,etpl_face,refine_flag); % Mesh refinement m. file

% Element polynomial smoothing algorithm
poly = smooth_func_poly(etpl_face,etpl.poly);
etpl.poly=poly;
