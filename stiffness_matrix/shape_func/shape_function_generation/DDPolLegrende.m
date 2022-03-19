function v = DDPolLegrende(x,p)
% Function to compute the second derivative of the Legendre polynomials
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

if p>=2 
   v(1:size(x,1),3) = 3; 
end

for i=3:p
    v(:,i+1) = ((2*i-1).*x.*v(:,i)-(1+i).*v(:,i-1))/(i-2);
end