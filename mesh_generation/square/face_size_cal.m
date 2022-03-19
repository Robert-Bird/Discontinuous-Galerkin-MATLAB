function [height]=face_size_cal(fn_loc,coord)
% Determines the size of an element face
%
% Input(s):
% fn_loc  - Local face number                                    
%
% coord   - a single elements coordinates 
%
% Ouput(s):
% height - lenght of face

%  Copyright (C) 2018 Robert Bird
%  $Revision: 1.0 $Date: 2018/06/11 17:09:20 $

face_node=[1 2;2 3;3 1];                                                   % Face to local node number conversion matrix
node_coord = coord(face_node(fn_loc,:),:);                                 % Positive element face coordinates
height = sqrt((node_coord(1,1)-node_coord(2,1))^2 ...                      % Element face length
            + (node_coord(1,2)-node_coord(2,2))^2);  
