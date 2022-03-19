function [tndof]=tndof_sum(element_p,nD)
% Determines the total number of degrees of freedom for in the volumetric integral
%
% Input(s):
% element_poly - Size of element polynomials 
% nD           - Number of dimensions
%
% Ouput(s):
% tndot        - Total number of degrees of freedom for volume integral

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

tndof=0;
for i=1:max(element_p(:,2))
    neDoF=(nov_calc(i)*nD)^2;
    tndof=tndof+neDoF*sum(element_p(:,2)==i);
end
