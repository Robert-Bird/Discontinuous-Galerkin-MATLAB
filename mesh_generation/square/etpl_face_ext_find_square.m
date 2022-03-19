function [etpl_face]=etpl_face_ext_find_square(coord,etpl,etpl_face,BC,conn)
% Defines the boundary coonditions on each face (Default homogeneous Neuman)
%
% Input(s):
% coord - Element coordinates                                         [x,y]
% etpl  - Element topology matrix                   [node 1,node 2, node 3]
% etpl_face  - Element face topology
% [element 1,element 2,local face 1,local face 2,normal x,normal y,face type]
%
% Ouput(s):
% etpl_face  - Element face topology
% [element 1,element 2,local face 1,local face 2,normal x,normal y,face type]
%              face type 1 - Internal element face
%
%                        2 - Homogeneous Dirichlet
%                       -2 - Prescribed  Dirichlet
%
%                        3 - Homogeneous Neumann
%                       -3 - Prescribed  Neumann
%
%                        4 - Homogeneous Dirichlet / Neumann
%                       -4 - Prescribed  Dirichlet / Neumann

%  Copyright (C) 2018 Robert Bird
%  $Revision: 1.0 $Date: 2018/06/11 17:09:20 $

etpl_face_int=etpl_face(etpl_face(:,2)~=0,:);                              % Internal faces
etpl_face_ext=etpl_face(etpl_face(:,2)==0,:);                              % External faces
fn_num=[1 2; 2 3; 3 1];                                                    % face to node numbering


for i = 1:size(etpl_face_ext,1)                                            % Searching through the external faces
    el_num=etpl_face_ext(i,1);
    loc_face_num=etpl_face_ext(i,3);
    n1=fn_num(loc_face_num,1);
    n2=fn_num(loc_face_num,2);
    node_1=etpl(el_num,n1);
    node_2=etpl(el_num,n2);
    
    index=(conn(:,1)==node_1).*(conn(:,2)==node_2)+(conn(:,2)==node_1).*(conn(:,1)==node_2);
    if sum(index)>0
        BC_expect=conn(index==1,3);
        etpl_face_ext(i,end)=conn(index==1,3);
        if BC_expect~=etpl_face_ext(i,end)
            fprintf('error in BC\n');
        end
    end
end

etpl_face=[etpl_face_int;etpl_face_ext];
