function [ngp,w,dNr,Nr]=dNr_vol_gauss_delta_global(el_num,poly_order,h,dx_flag,dy_flag,etpl,coord)
% Generates Shape functions, and their derivatives, for volumetric Gauss point locations
%
% Input(s):
% el_num     - Element number
% poly_order - polynomial order
% h          - space in h (positive or negative)
% dx_flag    - Flag to indicate change in x (==1)
% dy_flag    - Flag to indicate change in y (==1)
% etpl       - element topology structure
% coord      - list of coordinates in the mesh
%
% Ouput(s):
% w   - Gauss point weight values
% dNr - Shape function derivatives
% Nr  - Shape functions
% ngp - number of gauss points

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2018/06/11 17:09:20 $

pg=(4*(poly_order+1)); 
[xsi_s,eta_s,w]=Gauss_points(pg); 
[Nr]=local_sNr(xsi_s,eta_s,poly_order);
xsi_l=zeros(size(xsi_s));eta_l=xsi_s;
% Small element
for i = 1:size(xsi_s,1)
    el_coord_sm=[Nr(i,1:3)*coord(etpl.mat(el_num,:),:)]';                  % Finding the global coordinate
    glo_coord = [1;el_coord_sm];                                           % 
    
    if dx_flag==1
        glo_coord(2)=glo_coord(2)+h;                                       % If the dx_flag==1 add h to x
    elseif dy_flag==1
        glo_coord(3)=glo_coord(3)+h;                                       % If the dy_flag==1 add h to y
    end
    
    % Finding the adjusted Gauss point positions
    el_coord_lg=[1 1 1;coord(etpl.mat(el_num,:),:)'];
    Nx(1,:)=el_coord_lg\glo_coord;
    X_l=Nx*[-1 -1;1 -1;-1 1];
    xsi_l(i)=X_l(1);
    eta_l(i)=X_l(2);
end


[dNr]=local_dNr(xsi_l,eta_l,poly_order);                                   % Global Shape function derivatives
ngp=max(size(w));                                                          % Number of Gauss points


