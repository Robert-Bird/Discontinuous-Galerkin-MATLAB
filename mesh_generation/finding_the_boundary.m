function conn=finding_the_boundary(conn,BC,~,~)
% Function to match boundary conditions on edges.
%
% A tool, based on mesh topology, to match the boundary
% conditions of the input edges to refine2.m to the triangular mesh exterior
% output edges of refine2.m
%
% Input(s):
% conn       - Output edges of refine2
% BC         - Input edges to refine2 with the corresponding boundary condition
%
% Ouput(s):
% conn       - Updated edges

%  Copyright (C) 2018 Robert Bird
%  $Revision: 1.0 $Date: 2018/08/21 17:09:20 $

conn=[conn,3*ones(size(conn,1),1)];                                        % Expanding conn so the third row is the type of boundary condition (Default homogeneous Dirichlet)
seen=zeros(size(conn,1),1);                                                % Seen is an index saying if the element face (row in conn) has been visited already
BC_remaining=ones(size(BC,1),1);                                           % BC_remaining is an index of the BC conditions of the vertices have been assigned to the new mesh vertices    
count=0;                                                                   % Count keeps track of the boundary edges that have been assigned                                                        
original_nodes=unique(BC(:,1:2));                                          % Nodes numbers used as input into refine2 to create the new mesh
while count~=size(BC,1)                                                    % If count==size(BC,1) all the edges have been assigned
    exit = 0;i=0;                                                          % 'exit' - Exit flag and 'i' edge number of conn 
    while exit == 0                                                        % While loop to find an edge in conn that contains an original node
        i=i+1;                                                             % Go to the next edge
        found_index=conn(i,1)==original_nodes | conn(i,2)==original_nodes; % If the current edge 'i' contains an original node
        if sum(found_index)>0 && seen(i)==0                                % If the above line is true and edge 'i' hasn't been seen before
            index=[sum(conn(i,1)==original_nodes(found_index)),sum(conn(i,2)==original_nodes(found_index)),-1]; % Which nodes of the edge 'i' are a seed node
            start_node=conn(i,index==1);                                   % The first node of the remaining edges in conn that is a original node
            next_node=conn(i,index==0);                                    % The next node to look for in conn
            seen(i)=1;                                                     % Marking the edge 'i' as seen
            exit=1;                                                        % Exiting the loop
        end
    end
    
    if size(start_node)==1                                                 % If only one of the nodes in conn(i,1:2) is an original node
        faces_to_set=zeros(size(conn,1),1);                                % A vector to keep track of all the edges that are linked between two orginal node
        faces_to_set(i)=1;                                                 % Marking edge i as the first member in the chain
        exit = 0;ii=0;                                                     % Exit flag and counter through conn
        while exit == 0                                                    % Loop through all edges in conn until a chain of faces has been created between two original nodes
            ii=ii+1;                                                       % Increase counter
            if sum(conn(ii,1)==next_node | conn(ii,2)==next_node)>0 && seen(ii)==0 % If the edge ii contains next node, next edge in the chain is found
                index=[sum(conn(ii,1)==next_node),sum(conn(ii,2)==next_node),-1]; % Determine which node of edge i matches next_node
                next_node=conn(ii,index==0);                               % Determine the new next_node 
                faces_to_set(ii)=1;                                        % Add face ii to the chain
                seen(ii)=1;                                                % Mark edge ii as seen
            end
            if  sum(next_node==original_nodes)>0                           % If the next_node is an orginal node then the chain is complete
                exit=1;                                                    %
                end_node=next_node;                                        %
            end
            if ii == size(conn,1)                                          % If end of conn is reached loop start again
                ii=0;
            end
        end
        
        BC_index=(BC(:,1)== start_node &  BC(:,2)== end_node) | (BC(:,2)== start_node &  BC(:,1)== end_node); % Find the original face for the mesh generation which matches the start and end nodes
        BC_remaining(BC_index)=0;                                          % Make the original face as matched
        conn(faces_to_set==1,3)=BC(BC_index,3);                            % Set the boundary condition of all edges in faces_to_set to the same as the original edge that they lie on 
        count=count+1;                                                     
    else                                                                   % If generated edge matches an original edge then length(start)==2, and start = [start_node,end_node]
        BC_index=(BC(:,1)== start_node(1) &  BC(:,2)== start_node(2)) | (BC(:,2)== start_node(1) &  BC(:,1)== start_node(2)); 
        BC_remaining(BC_index)=0;
        conn(i,3)=BC(BC_index,3);
        count=count+1;
    end
end

end
