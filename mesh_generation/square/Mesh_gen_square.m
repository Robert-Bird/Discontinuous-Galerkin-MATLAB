function [etpl,coord,conn]=Mesh_gen_square(node,edge,BC)
% Mesh generation using triangle_mesh (a built-in MATLAB Delaunay-refinement
% mesher; a port of the relevant behaviour of Shewchuk's Triangle)
%
% Input(s):
% node     - Vertices of the desired domain (for mesh generation)
% edge     - Topology matrix of the desired domain (for mesh generation)
% BC       - Boundary conditions to be  imposed
% 
% Ouput(s):
% etpl  - Element topology matrix                   [node 1,node 2, node 3]
% coord - Element coordinates                                         [x,y]
% conn  - node pairs on the boundary and their corresponding BCs

%  Copyright (C) 2018 Robert Bird
%  $Revision: 1.0 $Date: 2018/06/11 17:09:20 $

hfun=0.1;                                                                  % Target element size (max edge length)
[vert,conn,tria] = triangle_mesh(node,edge,hfun);                          % Quality Delaunay mesh (built-in MATLAB port of Triangle)
conn = finding_the_boundary(conn,BC,tria,vert);
etpl=tria;
coord = vert;                                                              % Extracting nodal coordiantes data from text files
end







