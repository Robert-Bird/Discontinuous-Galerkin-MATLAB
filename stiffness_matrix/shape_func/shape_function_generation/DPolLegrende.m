function v = DPolLegrende(x,p)
% Function to compute the first derivative of the Legendre polynomials
%
% Input(s):
% x           - Evaluation points
% p           - Order
%
% Ouput(s):
% v           - Values of the derivatives

%  Copyright (C) 2017 Stefano Giani
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

v = zeros(size(x,1),p+1);

if p>=0 
   v(1:size(x,1),1) = 0; 
else
    assert(p>=0,'p must be non-negative');
end

if p>=1 
   v(1:size(x,1),2) = 1; 
end

for i=2:p
    v(:,i+1) = ((2*i-1).*x.*v(:,i)-i.*v(:,i-1))/(i-1);
end