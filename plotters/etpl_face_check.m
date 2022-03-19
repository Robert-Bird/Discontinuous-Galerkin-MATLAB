function etpl_face_check(etpl,coord,etpl_face,quiver_flag,face_numbers_flag,face_flag,element_numbers_flag)
% Plots the element topology and the element face topology
%
% Input(s):
% etpl                 - Element tolopogy struture
% coord                - Element coordinates
% etpl_face            - Element face tolopogy struture
% quiver_flag          - Plot the normal on each of the faces (1 = yes,0 = no)
% face_numbers_flag    - Flag to plot face numbers on mesh (1 = yes,0 = no)
% face_flag            - Plot element mesh (1 = yes,0 = no)
% element_numbers_flag - Flag to plot element numbers on mesh (1 = yes,0 = no)

%  Copyright (C) 2018 Thomas Wiltshire
%  $Revision: 1.0 $Date: 2018/08/21 17:09:20 $

figure                                                                     % Creating a new figure
index=zeros(size(etpl.mat,1),1);                                           % Index to check if an element is visited more than once for element numbering
hold on                                                                    %
face_to_nodes=[1 2;2 3;3 1];                                               % Face to nodes on face matrix
for i=1:size(etpl_face,1)                                                  % Looping through all faces in the mesh
    el=etpl_face(i,1);                                                     % Positive element for face i
    local_face=etpl_face(i,3);                                             % Local face of element el
    local_nodes=face_to_nodes(local_face,:);                               % Use face number to define local nodes
    global_nodes=etpl.mat(etpl_face(i,1),local_nodes);                     % Convert local nodes to global nodes
    X=coord(global_nodes,1);                                               % x,y coordinates of face nodes
    Y=coord(global_nodes,2);                                               %
    H = norm([X(2)-X(1);Y(2)-Y(1)]);                                       % Length of the current face
    
    if quiver_flag == 1                                                    %
        quiver(mean(X),mean(Y),etpl_face(i,5)*H/2,etpl_face(i,6)*H/2);     % Plotting the outward normal to the face i
    end                                                                    %
    
    if face_numbers_flag == 1                                              %
        s=sprintf('%.0d',i);                                               %
        text(mean(X),mean(Y),s);                                           % Printing global face number in plot
    end                                                                    %
    
    if face_flag == 1
        if    etpl_face(i,7)==1                                            % Boundary condition flag and plotting
            plot(X,Y,'k-')
        elseif    etpl_face(i,7)==2                                        % Boundary condition flag and plotting
            plot(X,Y,'r-')
        elseif etpl_face(i,7)==-2
            plot(X,Y,'b-')
        elseif etpl_face(i,7)==3
            plot(X,Y,'g-')
        elseif etpl_face(i,7)==-3
            plot(X,Y,'y-')
        elseif etpl_face(i,7)==4
            plot(X,Y,'m-')
        elseif etpl_face(i,7)==-4
            plot(X,Y,'c-')
        end
    end
    
    if element_numbers_flag==1 && index(el)==0                             % Plotting the element number
        index(el)=1;                                                       % Marking the element as seen so text print not repeated                                         
        s = sprintf('%.0f',i);
        text(mean(X),mean(Y),s);                                           % Printing global face number in plot
    end
end

if face_flag==1
    h(1) = plot(NaN,NaN,'k-');
    h(2) = plot(NaN,NaN,'r-');                                             % Plot nothing (for legend generation)
    h(3) = plot(NaN,NaN,'b-');
    h(4) = plot(NaN,NaN,'g-');
    h(5) = plot(NaN,NaN,'y-');
    h(6) = plot(NaN,NaN,'m-');
    h(7) = plot(NaN,NaN,'c-');
    legend(h,'Internal Face','Homogenous Dirichlet','Prescribed Dirichlet','Homogenous Neumann','Prescribed Neumann','Homogenous Dirichlet/ Neumann','Prescribed Dirichlet/Neumann');
end

axis equal                                                                 % Axis flags
axis off;