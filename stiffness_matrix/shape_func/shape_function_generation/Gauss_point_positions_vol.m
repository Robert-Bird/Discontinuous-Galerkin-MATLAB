function [Nr,dNr,vwp,ngp]=Gauss_point_positions_vol(p,el_sm,el_lg,coord,etpl)
% Maps the coordinates(xsi,eta) from: local(small element)->global->local(large element)
%
% Input(s):
% p       - Polynomial of the element
% el_sm   - Small element number
% el_lg   - Large element number
% coord   - coordinates of the entire mesh
% etpl    - Element topology structure
%
% Ouput(s):
% Nr - Shape functions
% dNr - Shape function derivatives
% vwp     - corresponding wieght of Gauss point
% ngp   - Number of Gauss points

%  Copyright (C) 2017 Robert Bird
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

pg=(4*(p+1));                                                              % Polynormal order or exact integration 1D
[xsi_s,eta_s,vwp]=Gauss_points(pg);                                              % Gauss points locations in the local domain and weights
dNr.s=local_dNr(xsi_s,eta_s,p);
Nr.s=local_sNr(xsi_s,eta_s,p);   

xsi_l=zeros(size(xsi_s));eta_l=zeros(size(xsi_s));

% Small element
for i = 1:size(vwp,1)
    el_coord_sm=[Nr.s(i,1:3)*coord(etpl(el_sm,:),:)]';
    glo_coord = [1;el_coord_sm];
    
    % large element
    el_coord_lg=[1 1 1;coord(etpl(el_lg,:),:)'];
    Nx(1,:)=el_coord_lg\glo_coord;
    X_l=Nx*[-1 -1;1 -1;-1 1];
    xsi_l(i)=X_l(1);
    eta_l(i)=X_l(2);
end

dNr.l=local_dNr(xsi_l,eta_l,p);Nr.l=local_sNr(xsi_l,eta_l,p); 
ngp=max(size(vwp));