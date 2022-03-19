function [pen,L]=pen_func(co_pve,co_nve,h,pos_loc_p,neg_loc_p,E,v,nx,ny)
% Creates the penalty value function for a given polynomial order
%
% Input(s):
% co_pve - Positive element coordinates
% co_nve - Negative element coordinates
% h      - Face length
% pos_loc_p - Positive element polynomial order
% neg_loc_n - Negative element polynomial order
% E      - Young's modulus
% v      - Poisson's ratio
% n_x    - x component of outward normal to face
% n_y    - y component of outward normal to face
%
% Ouput(s):
% pen   - Penalty value
% L     - Penalty term

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2018/06/11 17:09:20 $

a1=polyarea(co_pve(:,1),co_pve(:,2));
a2=polyarea(co_nve(:,1),co_nve(:,2));
a=min(a1,a2);                                                              % Extracting the smaller area
p=max(pos_loc_p,neg_loc_p);                                                % Maximum polynomial order
hmax=(2*a)/h;                                                              % Determining perpendicular length
N=[nx 0 ny;
   0 ny nx];
% D=N'*d_mat(E,v)*N;
D=d_mat(E,v);
D=max(eigs(D));

L=10*D;                                                                    % Penalty term
pen = ((p^2)*L)/hmax;                                                      % Complete coefficient
end
