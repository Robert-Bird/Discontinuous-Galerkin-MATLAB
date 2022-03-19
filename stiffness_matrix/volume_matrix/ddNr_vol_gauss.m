function [ngp,w,ddNr]=ddNr_vol_gauss(poly_order)
% Generates the second order derivatives of shape functions, for volumetric Gauss point locations
%
% Input(s):
% poly_order - Polynomial order
%
% Ouput(s):
% ngp  - Number of gauss points
% w    - Gauss point weight values
% ddNr - Second order derivatives in terms of xx and yy 

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

[xsi,eta,w]=Gauss_points(poly_order);                                              % Gauss points locations in the local domain and weights
[ddNr]=local_ddNr(xsi,eta,poly_order);                                     % Second order derivatives in terms of xx and yy  
ngp=max(size(w));                                                          % Number of Gauss points
