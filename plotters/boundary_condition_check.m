function boundary_condition_check(etpl,etpl_face,coord)
% Generates a colour coded plot of the boundary conditions applied to external faces. 
%
% Input(s):
% etpl       - Element topology structure
%             .mat:  Element topology matrix       
%             .poly: Element polynomial order list
%             .tree: Element tree structure
%
% etpl_face  - Element face topology                         
%
% coord      - Element coordinates

%  Copyright (C) 2018 Thomas Wiltshire
%  $Revision: 1.0 $Date: 2018/08/21 17:09:20 $

figure
triplot(etpl.mat,coord(:,1),coord(:,2),'k')                                % Create plot of mesh
face_to_nodes=[1 2;2 3;3 1];                                                
for i = 1:size(etpl_face,1)                                                % Step through each face in the mesh
    hold on;
    local_face=etpl_face(i,3);                                             % Fetch local face number
    local_nodes=face_to_nodes(local_face,:);                               % Use face number to define local nodes
    global_nodes=etpl.mat(etpl_face(i,1),local_nodes);                     % Convert local nodes to global nodes
    X=coord(global_nodes,1);                                               % x,y coordinates of face nodes
    Y=coord(global_nodes,2);
    if    etpl_face(i,7)==2                                                % Boundary condition flag and plotting
        plot(X,Y,'r-','LineWidth',2)
    elseif etpl_face(i,7)==-2
        plot(X,Y,'b-','LineWidth',2)
    elseif etpl_face(i,7)==3
        plot(X,Y,'g-','LineWidth',2)
    elseif etpl_face(i,7)==-3
        plot(X,Y,'y-','LineWidth',2)
    elseif etpl_face(i,7)==4
        plot(X,Y,'m-','LineWidth',2)
    elseif etpl_face(i,7)==-4
        plot(X,Y,'c-','LineWidth',2)
    end
end

hold on;
h = zeros(6, 1);
h(1) = plot(NaN,NaN,'r-');                                                 % Plot nothing (for legend generation)
h(2) = plot(NaN,NaN,'b-');
h(3) = plot(NaN,NaN,'g-');
h(4) = plot(NaN,NaN,'y-');
h(5) = plot(NaN,NaN,'m-');
h(6) = plot(NaN,NaN,'c-');

% Hard code legend colouring
legend(h, 'Homogenous Dirichlet','Prescribed Dirichlet','Homogenous Neumann','Prescribed Neumann','Homogenous Dirichlet/ Neumann','Prescribed Dirichlet/Neumann');

axis equal; axis off;
drawnow;