function [xsi,eta,w]=Gauss_points_s(p)
% Generates Gauss point locations for external element faces
%
% Input(s):
% p - polynomial order
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
xsi_f1 = (x+1).*(-1 -1)./2 + 1;
eta_f1 = (x+1).*(-1+1)./2 -1;

% face 3
xsi_f2 = (x+1).*(1+1)./2 -1;
eta_f2 = (x+1).*(-1 -1)./2 + 1;

% face 2
xsi_f3 = (x+1).*(-1+1)./2 - 1;
eta_f3 = (x+1).*(2)./2 + -1;

xsi=[xsi_f1;xsi_f2;xsi_f3];
eta=[eta_f1;eta_f2;eta_f3];
