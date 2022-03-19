function [tndof]=tndof_pp_f(element_p,nD,etpl_face)
% Determines the total number of degrees of freedom for surface integral face 1.
%
% Input(s):
% element_poly - Size of element polynomials 
% nD           - Number of dimensions
% etpl_face    - Element face tolopogy structure (see seed_mesh.m)
%
% Ouput(s):
% tnfot        - Total number of degrees of freedom for face 1

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

tndof=0;                                                                   % Settting tndof counter to zero.
for i=1:size(etpl_face,1)                                                  % Looping through the faces
    pe=element_p(etpl_face(i,1),2);                                        % Polynomial order of positive element
    p_nov=nov_calc(pe);                                                    % Calculates number of variables of positive elements
    neDoF=(p_nov*nD)^2;                                                    % Number of degrees of freedom for face i
    tndof=tndof+neDoF;                                                     % Increase tndof count
end
