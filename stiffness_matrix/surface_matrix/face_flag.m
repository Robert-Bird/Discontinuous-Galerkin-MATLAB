function [f_flag,h_small] = face_flag(etpl_face,coord_p,coord_n)
% Finds the smallest face between two elements
%
% Input(s):
% etpl_face     - Element topology of connected elements
% coord_p       - Coordinates of element etpl_face(1)
% coord_p       - Coordinates of element etpl_face(2)
%
% Ouput(s):
% f_flag        - Output if positive or negative face is smaller (used in debug)
% h_small       - Size of the smaller face

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2018/06/11 17:09:20 $

face_node = [1 2;2 3;3 1];
f(1)    = etpl_face(3); f(2)  = etpl_face(4);                          % local face number
if f(2)==0
    f(2)=f(1);
end
nds_1   = (face_node(f(1),:)); nds_2 = (face_node(f(2),:));
co_x_1  = coord_p(nds_1,1); co_x_2 = coord_n(nds_2,1);                     % face nodes
co_y_1  = coord_p(nds_1,2); co_y_2 = coord_n(nds_2,2);
h_1=norm([(co_x_1(1)-co_x_1(2)),(co_y_1(1)-co_y_1(2))]);                   % face length
h_2=norm([(co_x_2(1)-co_x_2(2)),(co_y_2(1)-co_y_2(2))]);
h_small=min([h_1,h_2]);
if h_1==h_2
    f_flag=[0,0];
else
    [~,f_flag]=sort([h_1,h_2]);
end
