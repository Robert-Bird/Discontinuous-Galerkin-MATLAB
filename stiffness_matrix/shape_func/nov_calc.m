function nov=nov_calc(poly_order)
% Returns the number of variables in 1D for the K triangle element.
% The method comes from: 
% P. Solin, K. Segeth, and I. Dolezel. Higher-order finite element methods. CRC Press, 2003.
%
% Input(s):
% poly_order - polynomial order
%
% Ouput(s):
% nov        - Number of variables in this element for 1D

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

if poly_order==2
    nov=12/2;
elseif poly_order>=3
    nov=(poly_order^2+3*poly_order+2)/2;
else
    nov=6/2;
end