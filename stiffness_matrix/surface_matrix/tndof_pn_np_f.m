function [tndof]=tndof_pn_np_f(element_poly,nD,etpl_face)
% Determines the total number of degrees of freedom for surface integral faces, 2 and 3.
%
% Input(s):
% element_poly - Size of element polynomials 
% nD           - Number of dimensions
% etpl_face    - Element face tolopogy structure (see seed_mesh.m)
%
% Ouput(s):
% tnfot        - Total number of degrees of freedom for faces 2 and 3

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

tndof=0;                                                                   % Settting tndof counter to zero.
for i=1:size(etpl_face,1)                                                  % Looping through the faces
    pe=element_poly(etpl_face(i,1),2);                                     % Polynomial order of positive element
    ne=element_poly(etpl_face(i,2),2);                                     % Polynomial order of negative element
    p_nov=nov_calc(pe);                                                    % Calculates number of variables of positive elements
    n_nov=nov_calc(ne);                                                    % Calculates number of variables of negative elements
    neDoF=p_nov*n_nov*(nD^2);                                              % Number of degrees of freedom for face i
    tndof=tndof+neDoF;                                                     % Increase tndof count
end