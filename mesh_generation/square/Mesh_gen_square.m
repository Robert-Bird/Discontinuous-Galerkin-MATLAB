function [etpl,coord,conn]=Mesh_gen_square(node,edge,BC)
% Mesh generation using Triangle c code
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

hfun=0.1;                                                                  % Mesh size function argument (see Mesh2D documentation)
[vert,conn,tria,~] = refine2(node,edge,[],[],hfun); 
conn = finding_the_boundary(conn,BC,tria,vert);
etpl=tria;
coord = vert;                                                              % Extracting nodal coordiantes data from text files
end







