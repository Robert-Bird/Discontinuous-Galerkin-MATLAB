function scaler_plot(etpl,coord,scaler_value)
% Plot a scalar variable.
%
% A surface plotter for the distribution of scaler values over
% a mesh. The scaler values is constant over an element such as the
% estimated error
%
% Input(s):
% etpl         - Element topology structure
% etpl_face    - Element face topology                         
% coord        - Element coordinates
% scaler_value - The scale value to plot

%  Copyright (C) 2018 Robert Bird
%  $Revision: 1.0 $Date: 2018/08/21 17:09:20 $

figure                                                                     % Creating figure
hold on                                                                    %
axis equal                                                                 % Axis flag
for i = 1:size(etpl.mat,1)                                                 % Looping over all elements in the mesh
   co=coord(etpl.mat(i,:),:);                                              % Nodes coordinates
   Er_n=repmat(scaler_value(i),3,1);                                       % Repeating the scaler value over all coordinates
   trisurf(1:3,co(:,1),co(:,2),Er_n);                                      % Plotting the element in 3D plot, z-axis is the scaler value 
end

