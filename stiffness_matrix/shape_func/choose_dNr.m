function dNr=choose_dNr(dNr,p,p_max)
% Return appropriate local shape function derivatives.
% Chooses local shape function derivatives given the maximum 
% plynomial order in the mesh and the polynomial order of the current element
%
% Input(s):
% dNr   - local shape function derivatives 
% p     - local polynomial order
% p_max - maximum polynomial order in the mesh 
%
% Ouput(s):
% dNr   -The reduced dNr matrix for the polynomial order p

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

% Selecting vertex shape functions
dNrv=dNr(:,1:3);

% Selecting surface shape functions
dNrs=dNr(:,(1:((p-1)*3))+3);

% Selecting bubble shape functions
b_start=((p_max-1)*3)+3;
dNrb=dNr(:,(b_start+1):end);
dNrb=dNrb(:,1:((p-1)*(p-2)/2));
dNr=[dNrv,dNrs,dNrb];
