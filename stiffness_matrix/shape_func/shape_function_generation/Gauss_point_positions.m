function [xsi_l,eta_l,xsi_s,eta_s,w,ngp]=Gauss_point_positions(p_max,el_sm,el_lg,coord,etpl,fp)
% Maps the coordinates(xsi,eta) from: local(small face)->global->local(large face)
%
% Input(s):
% p_max   - Maximum polynomial of elements sharing the face
% el_sm   - Small element number
% el_lg   - Large element number
% coord   - coordinates of the entire mesh
% etpl    - Element topology structure
% fp      - Local face number
%
% Ouput(s):
% xsi_l - Large element Gauss position xsi
% eta_l - Large element Gauss position eta 
% xsi_s - small element Gauss position xsi 
% eta_s - small element Gauss position eta 
% ngp   - Number of Gauss points
% w     - Gauss point weight 

%  Copyright (C) 2017 Robert Bird
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

[xsi_s,eta_s,w]=Gauss_points_s_mixed(p_max,fp);
xsi_l=zeros(size(xsi_s));eta_l=zeros(size(xsi_s));
[Nr.s]=local_sNr(xsi_s,eta_s,p_max);

% Small element
for i = 1:size(xsi_s,1)
    el_coord_sm=[Nr.s(i,1:3)*coord(etpl(el_sm,:),:)]';
    glo_coord = [1;el_coord_sm];
    
    % large element
    el_coord_lg=[1 1 1;coord(etpl(el_lg,:),:)'];
    Nx(1,:)=el_coord_lg\glo_coord;
    X_l=Nx*[-1 -1;1 -1;-1 1];
    xsi_l(i)=X_l(1);
    eta_l(i)=X_l(2);
end

ngp=max(size(w));