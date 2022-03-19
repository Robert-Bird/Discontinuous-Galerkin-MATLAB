function ddNr=choose_ddNr_xy(ddNr,p,p_max)
% Return appropriate local shape function second derivatives.
% Chooses local shape function (2nd order) derivatives given 
% the maximum  polynomial order in the mesh and the polynomial order of 
% the current element
%
% Input(s):
% ddNr   - local shape function second derivatives 
% p      - local polynomial order
% p_max  - maximum polynomial order in the mesh 
%
% Ouput(s):
% ddNr   -The reduced ddNr matrix for the polynomial order p

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

% Selecting vertex shape functions
ddNrv=ddNr(:,1:3);

% Selecting surface shape functions
ddNrs=ddNr(:,(1:((p-1)*3))+3);

% Selecting bubble shape functions
b_start=((p_max-1)*3)+3;
ddNrb=ddNr(:,(b_start+1):end);
ddNrb=ddNrb(:,1:((p-1)*(p-2)/2));
ddNr=[ddNrv,ddNrs,ddNrb];
