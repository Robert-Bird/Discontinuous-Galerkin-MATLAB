function [coord,etpl,etpl_face]=seed_mesh_square(node,edge,BC)
% Main function for generating DG mesh topology
%
% Input(s):
% node     - Vertices of the desired domain (for mesh generation)
% edge     - Topology matrix of the desired domain (for mesh generation)
% BC       - Boundary conditions to be  imposed
%
% Ouput(s):
% coord      - Element coordinates                                    [x,y]
% 
% etpl       - Element topology structure
%             .mat:  Element topology matrix        [node 1,node 2, node 3]
%
%             .poly: Element polynomial order list
%                                         [element number,polynomial order]
%
%             .tree: Element tree structure
%                          [active(yes=1),generation number,parent element]
%
% etpl_face  - Element face topology

%  Copyright (C) 2018 Robert Bird
%  $Revision: 1.0 $Date: 2018/06/11 17:09:20 $

[etpl_mat,coord,conn]=Mesh_gen_square(node,edge,BC);                       % Square unstructured mesh
etpl.mat=etpl_mat;                                                         %
etpl.poly=[[1:size(etpl.mat,1)]',3*ones(size(etpl.mat,1),1)];              % Polynomial order
etpl_face=etpl_face_square_func(coord,etpl,BC,conn);                       % Generate etpl_face
etpl.tree=[ones(size(etpl.mat,1),1)==1,zeros(size(etpl.mat,1),1),zeros(length(etpl.mat(:,1)),1)]; % Element tree






