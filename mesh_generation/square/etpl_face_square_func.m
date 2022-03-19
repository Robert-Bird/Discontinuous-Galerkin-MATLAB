function  etpl_face=etpl_face_square_func(coord,etpl,BC,conn)
% Extracts element face topology from etpl and coord
%
% Input(s):
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
%             .ed_recale: List of elements to recalculate local stiffness
%                                                                [1=recalc]
%
% Ouput(s):
% etpl_face  - Element face topology
% [element 1,element 2,local face 1,local face 2,normal x,normal y,face type]

%  Copyright (C) 2018 Robert Bird
%  $Revision: 1.0 $Date: 2018/06/11 17:09:20 $

etpl_face = etpl_face_all(etpl,coord);                                     % All face information
etpl_face = etpl_face_ext_find_square(coord,etpl.mat,etpl_face,BC,conn);             % Added external face information
