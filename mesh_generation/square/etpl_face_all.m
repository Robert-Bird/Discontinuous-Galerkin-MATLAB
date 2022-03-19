function etpl_face = etpl_face_all(etpl,coord)
% Determines all face topology
%
% Input(s):
% etpl       - Element topology structure
%             .mat:  Element topology matrix        [node 1,node 2, node 3]
%
%             .poly: Element polynomial order list
%                                         [element number,polynomial order]
%
%             .tree: Element tree structure
%                          [active(yes=1),generation number,parent element]
%
% coord      - Element coordinates
%
% Ouput(s):
% etpl_face  - Element face topology
% [element 1,element 2,local face 1,local face 2,normal x,normal y,face type]

%  Copyright (C) 2018 Robert Bird
%  $Revision: 1.0 $Date: 2018/06/11 17:09:20 $

etpl_face_index=ones(size(etpl.mat));                                      % Creating index to keep track of faces already visited
face_index=[1 2;2 3;3 1];                                                  % Node index for faces: 1,2 and 3
etpl_face=zeros(size(etpl.mat,1)*3,7);                                     % Allocating memory for etpl_face
face_count=0;                                                              % Row counter for etpl_face
el_list(:,1)=1:size(etpl.mat,1);                                           % List of all elements in the mesh
face_list=1:3;                                                             
etpl_face_nodes=[etpl.mat(:,[1 2]),etpl.mat(:,[2 3]),etpl.mat(:,[3 1])];   % Matrix of nodes for faces 1,2 and 3 anticlockwise

for el = 1:size(etpl_face_index,1)                                         % Looping through all elements
    for loc_face = 1:size(etpl_face_index,2)                               % Looping through element faces
        if etpl_face_index(el,loc_face)==1                                 % If the face has not been used before
            pos_nodes=etpl_face_nodes(el,(loc_face*2-1:loc_face*2));       % Positive element nodes
            face_count=face_count+1;                                       % Face topology row count
            etpl_face(face_count,[1,3,7])=[el,loc_face,3];                 % Adding positive element characteristics to the matrix
            
            p_co = coord(etpl.mat(el,:),:);                                % Positive element nodes
            pf_c = p_co(face_index(loc_face,:),:);                         % Positive element face coordinates
            h = sqrt((pf_c(1,1)-pf_c(2,1))^2 + (pf_c(1,2)-pf_c(2,2))^2);   % Length of face
            n= [0 -1;1 0]*[(pf_c(1,1)-pf_c(2,1));pf_c(1,2)-pf_c(2,2)]./h;  % Outward normal to the face
            etpl_face(face_count,[5,6])=n;                                 % Adding normal to element face topology matrix
            
            el_neg   = el_list(sum(etpl_face_nodes(:,2:2:end)==pos_nodes(1),2)==1); % Finding possible negative element
            for i = 1:length(el_neg)                                       % Looping through possible negative elements
                neg_face=face_list(etpl_face_nodes(el_neg(i),2:2:end)==pos_nodes(1)); % Negative element face
                neg_global_node=etpl_face_nodes(el_neg(i),(neg_face*2)-1); % Corressponding node to negative element face
                if neg_global_node==pos_nodes(2)                           % Confirming that the positive and element elements have an adjacent face
                    etpl_face(face_count,[2,4,7])=[el_neg(i),neg_face,1];  % Adding negative element contributions to etpl_face
                    etpl_face_index(el_neg(i),neg_face)=0;                 % Face is marked as visited
               end
            end
        end
    end
end

% Reorganising etpl_face, lists internal first then external
internal=etpl_face(etpl_face(:,2)~=0,:);
external=etpl_face(etpl_face(:,2)==0,:);
etpl_face=[internal;external];
etpl_face(etpl_face(:,1)==0,:)=[];
