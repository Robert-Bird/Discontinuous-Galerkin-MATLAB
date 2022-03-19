function [tndof]=tndof_nn_f(element_p,nD,etpl_face)
% Determines the total number of degrees of freedom for surface integral face 4.
%
% Input(s):
% element_poly - Size of element polynomials 
% nD           - Number of dimensions
% etpl_face    - Element face tolopogy structure (see seed_mesh.m)
%
% Ouput(s):
% tnfot        - Total number of degrees of freedom for face 4

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

tndof=0;                                                                   % Settting tndof counter to zero.
for i=1:size(etpl_face,1)                                                  % Looping through the faces
    ne=element_p(etpl_face(i,2),2);                                        % Polynomial order of positive element
    n_nov=nov_calc(ne);                                                    % Calculates number of variables of positive elements
    neDoF=(n_nov*nD)^2;                                                    % Number of degrees of freedom for face i
    tndof=tndof+neDoF;                                                     % Increase tndof count
end
