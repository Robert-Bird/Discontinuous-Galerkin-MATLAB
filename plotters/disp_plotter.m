function disp_plotter(etpl,coord,ed,u,s)
% Creates nodal displacement plot of the mesh 
%
% Input(s):
% etpl       - Element tolopogy struture 
% coord      - Element coordinates 
% ed         - Degrees of freedom steering matrix
% u          - Displacement solution created by the SIPG algorithm
% s          - Scale factor

%  Copyright (C) 2018 Robert Bird
%  $Revision: 1.0 $Date: 2018/08/21 17:09:20 $

figure;                                                                    % Create a new figure
hold on;                                                                   % keep all lines plotted
for i = 1: size(etpl.mat,1)                                                % Looping through all elements
    if etpl.tree(i,1)==1                                                   % If element is in current mesh
        plot([coord(etpl.mat(i,:),1)+u(ed(i,1:2:6))*s;coord(etpl.mat(i,1),1)+u(ed(i,1))*s],...
            [coord(etpl.mat(i,:),2)+u(ed(i,2:2:6))*s;coord(etpl.mat(i,1),2)+u(ed(i,2))*s],...
            'k');                                                          % Plot a triangle's vertex displacement as three lines
    end
end
axis equal                                                                 % Axis flags
axis off                                                                   %