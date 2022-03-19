function [w,dNr,Nr,ngp]=dNr_surf_gauss(poly_order)
% Generates Shape functions, and their derivatives, for surface Gauss point locations on external faces
%
% Input(s):
% poly_order - polynomial order
%
% Ouput(s):
% w   - Gauss point weight values
% dNr - Shape function derivatives
% Nr  - Shape functions
% ngp - number of gauss points

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2018/06/11 17:09:20 $

[xsi,eta,w]=Gauss_points_s(poly_order);                              % Gauss point locations
[dNr]=local_dNr(xsi,eta,poly_order);                                       % Shape function derivatives
[Nr]=local_sNr(xsi,eta,poly_order);                                        % Shape functions
ngp=max(size(w));                                                          % number of gauss points


