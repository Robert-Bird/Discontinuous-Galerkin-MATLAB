function [ngp,w,dNr,Nr]=dNr_vol_gauss(poly_order)
% Generates Shape functions, and their derivatives, for volumetric Gauss point locations
%
% Input(s):
% poly_order - Polynomial order
%
% Ouput(s):
% ngp - number of gauss points
% w   - Gauss point weight values
% dNr - Shape function derivatives
% Nr  - Shape functions

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

[xsi,eta,w]=Gauss_points(poly_order);                                              % Gauss points locations in the local domain and weights
[dNr]=local_dNr(xsi,eta,poly_order);[Nr]=local_sNr(xsi,eta,poly_order);    % Global shape functions
ngp=max(size(w));                                                          % Number of Gauss points


