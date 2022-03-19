function [el_p2,el_n2,f_p2,f_n2,f]=pos_neg_el_switch(coord_p,coord_n,el_p,el_n,f_p,f_n)
% Switches faces if negative element has a smaller face
%
% Input(s):
% coord_p - Positive element coordinates
% coord_n - Negative element coordinates
% el_p    - Positive element
% el_n    - Negative element
% f_p     - Positive element face number
% f_n     - Negative element face number
%
% Ouput(s):
% el_p2   - New positive element
% el_n2   - New negative element
% f_p2    - New Positive element face number
% f_n2    - New Negative element face number
% f       - Original orientation of faces numbers [positive, negative]

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

face_node=[1 2;2 3;3 1];                                                   % Face to local node number conversion matrix
f(1)    = f_p; f(2)  = f_n;                                                % Local face numbers
nds_1   = (face_node(f(1),:)); nds_2 = (face_node(f(2),:));                % Local nodes of face numbers
co_x_1  = coord_p(nds_1,1); co_x_2 = coord_n(nds_2,1);                     % x coordinates of nodes on faces
co_y_1  = coord_p(nds_1,2); co_y_2 = coord_n(nds_2,2);                     % x coordinates of nodes on faces
h_1=norm([(co_x_1(1)-co_x_1(2)),(co_y_1(1)-co_y_1(2))]);                   % face length for positive element
h_2=norm([(co_x_2(1)-co_x_2(2)),(co_y_2(1)-co_y_2(2))]);                   % face length for negative element
[~,index]=sort([h_1,h_2]);                                                 % Find which face is smaller
if index(1)==1                                                             % Do nothing
    f=1;
    el_p2=el_p;
    el_n2=el_n;
    f_p2 =f_p;
    f_n2 =f_n;
else                                                                       % reverse faces
    f=-1;
    el_p2=el_n;
    el_n2=el_p;
    f_p2 =f_n;
    f_n2 =f_p;
end
