function De=d_mat(E,v)
% Function to compute the stiffness tensor
%
% Input(s):
% E          - Young's modulus
% v          - Poisson's ratio
%
% Ouput(s):
% De     - Stiffness tensor

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $
De=E/((1+v)*(1-(2*v)))*[(1-v) v 0;v (1-v) 0;0 0 (1-(2*v))/2];