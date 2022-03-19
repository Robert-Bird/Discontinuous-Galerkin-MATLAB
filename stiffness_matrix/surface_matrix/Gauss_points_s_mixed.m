function [xsi,eta,w]=Gauss_points_s_mixed(p,fp)
% Generates Gauss point locations for internal element faces
%
% Input(s):
% poly_order - polynomial order
% face_num   - face number
%
% Ouput(s):
% xsi - xsi locations of Gauss point for a face
% eta - eta locations of Gauss point for a face
% w   - corresponding wieght of Gauss point

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2018/06/11 17:09:20 $

n=(p+1);
q=(1:n)./(1:2:2*n);
e=(1:n-1)./(3:2:2*n);
b0=2;
[~,D,V]=svd(diag(sqrt(q))+diag(sqrt(e),1));
x=(diag(D).^2)-1;
w=b0*V(1,:).^2;

% face 1
if fp==1
    xsi = (x+1).*(-1-1)./2 + 1;
    eta = (x+1).*(-1+1)./2 -1;
elseif fp==2
    xsi = (x+1).*(1+1)./2 -1;
    eta = (x+1).*(-1-1)./2 + 1;
elseif fp==3
    xsi = (x+1).*(-1+1)./2 - 1;
    eta = (x+1).*(2)./2 + -1;
end
