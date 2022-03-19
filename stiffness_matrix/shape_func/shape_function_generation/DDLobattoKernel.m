function v = DDLobattoKernel(x,p)
% Function to compute the first second of the Lobatto kernel functions
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

assert(p>=0,'p must be non-negative');


l = DDDPolLegrende(x,p+1);


for i=2:p+2
    v(:,i-1) = -4*(2*i-1)*sqrt(i-0.5)/(i*(i-1))*l(:,i);
end