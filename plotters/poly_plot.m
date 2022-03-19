function poly_plot(etpl,coord)
% Plots a grey-scale colormap of polynomial order distribution in the mesh
%
% Input(s):
% etpl       - Element tolopogy struture 
% coord      - Element coordinates 

%  Copyright (C) 2018 Robert Bird
%  $Revision: 1.0 $Date: 2018/08/21 17:09:20 $

figure;                                                                    % Creates a figure                                    
hold on;                                                                   % 
for i=1:size(etpl.mat,1)                                                   % Loop through all elements used
    if etpl.tree(i,1)==1                                                   % Check if element is in the current mesh
        nodes=etpl.mat(i,:);                                               % nodes of element i
        x=coord(nodes,1);                                                  % x coordinate of element vertices
        y=coord(nodes,2);                                                  % y coordinate of element vertices
        z=repmat(etpl.poly(i,2),3,1);                                      % Create a z coordinate corresponding to polynomial order of element i
        trisurf(1:3,x,y,z);                                                % Plotting element with associated polynomial order
    end
end
axis equal                                                                 % Axis flags
axis off                                                                   %

                                                                           % Making the image a discrete greyscale colourmap
min_poly=min(etpl.poly(etpl.tree(:,1)==1,2));                              % Minimum polynomial order
max_poly=max(etpl.poly(etpl.tree(:,1)==1,2));                              % Minimum polynomial order
poly_range = max_poly-min_poly+1;                                          % 
gray_colour=gray;                                                          %
jump_gray=ceil(size(gray,1)./poly_range);                                  %
gray_colour=gray_colour(end:-jump_gray:1,:);                               % Making the number of colours in the colour map to the range in polynomial order
colorbar                                                                   %
colormap(gray_colour);                                                     % Changing the colourmap to the new discrete colourmap

ax=gca;
uisetfont(ax);