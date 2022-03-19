function [coord,Nn4,Nn5,Nn6,max_node]=coord_check_new(nel,etpl,etpl_face,coord,N4,N5,N6,max_node)
% Function to check if a new node exists or not in the refinement tree.
%
% Input(s):
% nel         - Current element
% etpl        - Element topology structure
% etpl_face   - Element face topology matrix
% coord       - Element coordinates
% N4          - Mid-point node number between 1 and 2
% N5          - Mid-point node number between 2 and 3 
% N6          - Mid-point node number between 3 and 1 
% max_node    - Mumber of nodes in the mesh
%
% Ouput(s):
% coord       - Updated element coordinates
% Nn4          - Updated mid-point node number between 1 and 2
% Nn5          - Updated mid-point node number between 2 and 3 
% Nn6          - Updated mid-point node number between 3 and 1 
% Max_node    - Updated number of nodes in the mesh

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

coordt=coord(1:max_node,:);                                                % Temporary storing all the coodinates in the mesh
new_nodes=zeros(3,2);                                                      % Creating a matric of new nodes
new_node_count=0;                                                          % List of node numbers
face_list=[2 3 1];                                                         % Refined element 4 parent face -> node numbers
el_list=1:size(etpl.mat);                                                  % List of elements in the mesh

% Rearrange etpl_face into a more convientient form for faces connected to
% current parent element 'nel' in inputs, such that etpl_face_temp has the
% form:
%  parent face number | opposite element number |  opposite element face
%  number
nel_as_positive_opposite_element=etpl_face(etpl_face(:,1)==nel,2);
nel_as_positive_opposite_face=etpl_face(etpl_face(:,1)==nel,4);
nel_as_positive_face=etpl_face(etpl_face(:,1)==nel,3);
etpl_face_temp=[nel_as_positive_face,nel_as_positive_opposite_element,nel_as_positive_opposite_face];

nel_as_negative_opposite_element=etpl_face(etpl_face(:,2)==nel,1);
nel_as_negative_opposite_face=etpl_face(etpl_face(:,2)==nel,3);
nel_as_negative_face=etpl_face(etpl_face(:,2)==nel,4);
etpl_face_temp=[etpl_face_temp;[nel_as_negative_face,nel_as_negative_opposite_element,nel_as_negative_opposite_face]];

% Face 1 ------------------------------------------------------------------
face_1_info=etpl_face_temp(etpl_face_temp(:,1)==1,:);
if size(face_1_info,1)>1                                                   % Check the nodal coordinates and whether they exist already
    parent=etpl.tree(face_1_info(1,2),3);
    el_4=max(el_list(etpl.tree(:,3)==parent));
    common_face=face_1_info(1,3);
    Nn4=etpl.mat(el_4,face_list(common_face));                             % Find the already existing node
else
    max_node=max_node+1;                                                   % If the node does not exist then create a new node
    new_node_count=new_node_count+1;
    new_nodes(new_node_count,:)=N4;
    Nn4=max_node;
end

% Face 2 ------------------------------------------------------------------
face_2_info=etpl_face_temp(etpl_face_temp(:,1)==2,:);
if size(face_2_info,1)>1                                                   % Check the nodal coordinates and whether they exist already
    parent=etpl.tree(face_2_info(1,2),3);
    el_4=max(el_list(etpl.tree(:,3)==parent));
    common_face=face_2_info(1,3);
    Nn5=etpl.mat(el_4,face_list(common_face));                             % Find the already existing node
else
    max_node=max_node+1;                                                   % If the node does not exist then create a new node
    new_node_count=new_node_count+1;
    new_nodes(new_node_count,:)=N5;
    Nn5=max_node;
end

% Face 3 ------------------------------------------------------------------
face_3_info=etpl_face_temp(etpl_face_temp(:,1)==3,:);
if size(face_3_info,1)>1                                                   % Check the nodal coordinates and whether they exist already
    parent=etpl.tree(face_3_info(1,2),3);
    el_4=max(el_list(etpl.tree(:,3)==parent));
    common_face=face_3_info(1,3);
    Nn6=etpl.mat(el_4,face_list(common_face));                             % Find the already existing node
else
    max_node=max_node+1;                                                   % If the node does not exist then create a new node
    new_node_count=new_node_count+1;
    new_nodes(new_node_count,:)=N6;
    Nn6=max_node;
end
coord(1:max_node,:)=[coordt;new_nodes(1:new_node_count,:)];                % Update the coordinates