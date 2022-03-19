function [coord,etpl,etpl_face]=mesh_refine_topology_marking(coord,etpl,etpl_face,h_refine)
% Function to adapt the mesh in h and change the corresponding data structure.
%
% Input(s):
% coord       - Element coordinates
% etpl        - Element topology structure
% etpl_face   - Element face topology matrix
% h_refine    - Elements to be refined in h
%
% Ouput(s):
% coord       - Updated element coordinates
% etpl        - Updated element topology structure
% etpl_face   - Updated element face topology matrix
%
%  See also H_ADAPT.

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

% Setting up variables for changes to the element topology data structure
tree_counter=size(etpl.tree,1);                                            % Total number of active and unactive elements
fn_loc = [2,3,1];                                                          % Face number for interior element against siblings 1 2 3
local_face_list=[1 2;2 3;3 1];                                             % Nodes for faces 1,2 and 3
max_node=max(etpl.mat(:));                                                 % Current number of nodes in the mesh
nel=size(etpl.mat,1);                                                      % Current number of elements in the mesh

etpl.tree=[etpl.tree;zeros(tree_counter*4,3)];                             % Allocating space for new rows of etpl.tree
etpl.mat =[etpl.mat; zeros(size(etpl.mat,1)*4,3)];                         % Allocating space for new element topology
etpl.poly=[etpl.poly;zeros(size(etpl.poly,1)*4,2)];                        % Allocating space for new element polynomial order
coord=[coord;zeros(size(coord,1)*3,2)];                                    % Allocating space for new coordinates

etpl_face_counter=size(etpl_face,1);                                       % Counter for etpl_face
etpl_face=[etpl_face;zeros(size(etpl_face,1)*4,7)];                        % Allocating space for etpl_face

elements_to_refine=1:nel;                                                  % List of elements to refine
elements_to_refine=elements_to_refine(h_refine(:,2)==1);

h_refine_temp=elements_to_refine;
for i=0:max(etpl.tree(h_refine_temp,2))
    for k=1:length(h_refine_temp)
        if etpl.tree(h_refine_temp(k),2)==i
            % refine etpl
            ns = etpl.mat(h_refine_temp(k),:);% First select the nodes from etpl
            
            % Nodal coordinates
            N1=[coord(ns(1),1),coord(ns(1),2)];                        % Node 1 of parent element
            N2=[coord(ns(2),1),coord(ns(2),2)];                        % Node 2 of parent element
            N3=[coord(ns(3),1),coord(ns(3),2)];                        % Node 3 of parent element
            
            % New nodal coordinates, middle of edges
            N4=(N1+N2)./2;                                             % Node 4, midpoint between nodes 1 and 2
            N5=(N2+N3)./2;                                             % Node 5, midpoint between nodes 2 and 3
            N6=(N3+N1)./2;                                             % Node 6, midpoint between nodes 3 and 1
            
            % node numbers
            Nn1=ns(1);                                                 % Global node numbers
            Nn2=ns(2);                                                 % Global node numbers
            Nn3=ns(3);                                                 % Global node numbers
            
            % Determining if nodes N4, N5 N6 exist of not
            [coord,Nn4,Nn5,Nn6,max_node]=coord_check_new(h_refine_temp(k),etpl,etpl_face(etpl_face(:,1)~=0,:),coord,N4,N5,N6,max_node);
            
            % Creating new elements
            etpl.mat(nel+1,:)=[Nn1 Nn4 Nn6]; % face 1 and 3 exterior
            etpl.mat(nel+2,:)=[Nn4 Nn2 Nn5]; % face 1 and 2 exterior
            etpl.mat(nel+3,:)=[Nn6 Nn5 Nn3]; % face 2 and 3 exterior
            etpl.mat(nel+4,:)=[Nn6 Nn4 Nn5]; % Interior element
            
            % Updating polynomial order of elements
            etpl.poly((1:4)+nel,:)=[(1:4)+nel;repmat(etpl.poly(h_refine_temp(k),2),1,4)]'; % Child elements are same polynomial order as parents
            
            % Updating element tree
            etpl.tree(tree_counter+(1:4),:)=repmat([1,etpl.tree(h_refine_temp(k),2)+1,h_refine_temp(k)],4,1); % Creating new rows in the tree for children
            etpl.tree(h_refine_temp(k),1)=0;                           % Turning off parent element
            tree_counter=tree_counter+4;                               % Increasing the tree counter
            
            
            exterior=nel+(1:3);                                        % Determining the exterior children elements
            interior=nel+4;                                            % Determining the interior child element
            nel=nel+4;
            
            
            % Extracting etpl_face for the current parent element and defining as positive to all its adjacent elements for convienince
            etpl_face_current_el=[etpl_face(etpl_face(:,1)==h_refine_temp(k),:);etpl_face(etpl_face(:,2)==h_refine_temp(k),[2 1 4 3]),-etpl_face(etpl_face(:,2)==h_refine_temp(k),5:6),etpl_face(etpl_face(:,2)==h_refine_temp(k),7)];
            
            % Removing etpl_face corressponding to the parent element from etpl_face
            etpl_face(etpl_face(:,1)==h_refine_temp(k),1:2)=zeros(sum(etpl_face(:,1)==h_refine_temp(k),1),2);
            etpl_face(etpl_face(:,2)==h_refine_temp(k),1:2)=zeros(sum(etpl_face(:,2)==h_refine_temp(k),1),2);
            
            % Looping throught the three faces of the parent element
            for j=1:3
                etpl_face_current_face = etpl_face_current_el(etpl_face_current_el(:,3)==j,:);
                for t = size(etpl_face_current_face,1)
                    if size(etpl_face_current_face,1) == 1 % Both elements are initially the same size or an exterior
                        new_els = exterior(local_face_list(j,:));
                        etpl_face(etpl_face_counter+1,:)=[new_els(1),etpl_face_current_face(2:end)];
                        etpl_face(etpl_face_counter+2,:)=[new_els(2),etpl_face_current_face(2:end)];
                        etpl_face_counter=etpl_face_counter+2;
                    elseif  size(etpl_face_current_face,1) == 2 % The parent element is connect to two adjacent elements
                        for i1 =1:2
                            el_1=etpl_face_current_face(i1,2);             % Adjacent element selection
                            face_1=etpl_face_current_face(i1,4);           % Adjacent face selection
                            nodes_negative=etpl.mat(el_1,local_face_list(face_1,[2,1])); % Adjacent element nodes (but reversed for the current face)
                            nodes_positive=etpl.mat(h_refine_temp(k),local_face_list(j,:)); %  parent element nodes
                            
                            % As parent shares only 1 node with the
                            % adjacent smaller element (el_1), need to determine
                            % the node of the parent that is connected to
                            % el_1 and therefore the correect child of the
                            % parent that is connected to el_1
                            if nodes_positive(1)==nodes_negative(1)         
                                etpl_face_counter=etpl_face_counter+1;
                                
                                % Adding row for element face topology for the child element 
                                etpl_face(etpl_face_counter,:)=[exterior(local_face_list(j,1)),etpl_face_current_face(i1,2:end)];
                            elseif nodes_positive(2)==nodes_negative(2)
                                etpl_face_counter=etpl_face_counter+1;
                                
                                % Adding row for element face topology for the child element 
                                etpl_face(etpl_face_counter,:)=[exterior(local_face_list(j,2)),etpl_face_current_face(i1,2:end)];
                            end
                            
                        end
                    elseif  size(etpl_face_current_face,1) > 2
                        fprintf('Error in mesh refinement check mesh_refine_topology_marking.m\n')
                    end
                end
            end
            
            % Include the interior element to the new refinement
            for j2 = 1:3
                etpl_face_counter=etpl_face_counter+1;
                co_pve  = coord(etpl.mat(interior,:),:);                               % Positive coordinates
                pf_co = co_pve(local_face_list(j2,:),:);                               % Positive element face coordinates
                h = sqrt((pf_co(1,1)-pf_co(2,1))^2 + (pf_co(1,2)-pf_co(2,2))^2);       % Element face length
                n= [0 -1;1 0]*[(pf_co(1,1)-pf_co(2,1));pf_co(1,2)-pf_co(2,2)]./h;      % Outward normal
                etpl_face(etpl_face_counter,:)=[interior,exterior(j2),j2,fn_loc(j2),n(1),n(2),1];
            end
            
        end
    end
end

% Removing excess rows or deleted rows that are not required
etpl.tree((tree_counter+1):end,:)=[];
coord((max_node+1):end,:)=[];
etpl_face(etpl_face(:,1)==0,:)=[];
etpl.poly(etpl.mat(:,1)==0,:)=[];
etpl.mat(etpl.mat(:,1)==0,:)=[];
end










