function [dNr,Nr]=dNr_surf_gauss_mixed(max_poly,xsi_s,eta_s,xsi_l,eta_l)
% Determines the shape function values for the small and large faces
%
% Input(s):
% max_poly - Maximum polynomial order of elements sharing the face
% xsi_s    - xsi values for the small face
% eta_s    - eta values for the small face
% xsi_l    - xsi values for the alrge face
% eta_l    - eta values for the large face
%
% Ouput(s):
% dNr   - Shape function derivative structure
%       .p small face 
%       .n small face 
% Nr    - Shape function structure
%       .p small face 
%       .n small face 

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2018/06/11 17:09:20 $

[Nr.s]=local_sNr(xsi_s,eta_s,max_poly);                                    % Shape functions for small element
[dNr.s]=local_dNr(xsi_s,eta_s,max_poly);                                   % Shape functions derivatives for small element
[Nr.l]=local_sNr(xsi_l,eta_l,max_poly);                                    % Shape functions for large element
[dNr.l]=local_dNr(xsi_l,eta_l,max_poly);                                   % Shape functions derivatives for large element
dNr.p=dNr.s;                                                               % Small element is postive     
Nr.p=Nr.s;                                                                 %
dNr.n=dNr.l;                                                               % large element is negative    
Nr.n=Nr.l;                                                                 %   
clear dNr.l dNr.s Nr.l Nr.s                                                % Clearing unecessary variables


