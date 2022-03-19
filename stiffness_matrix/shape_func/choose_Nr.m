function Nr=choose_Nr(Nr,p,p_max)
% Return appropriate local shape function.
% Chooses local shape functions given the maximum 
% polynomial order in the mesh and the polynomial order of the current element
%
% Input(s):
% Nr    - local shape function derivatives 
% p     - local polynomial order
% p_max - maximum polynomial order in the mesh 
%
% Ouput(s):
% Nr   -The reduced Nr matrix for the polynomial order p

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

% Selecting vertex shape functions
Nrv=Nr(:,1:3);

% Selecting surface shape functions
Nrs=Nr(:,(1:((p-1)*3))+3);

% Selecting bubble shape functions
b_start=((p_max-1)*3)+3;
Nrb=Nr(:,(b_start+1):end);
Nrb=Nrb(:,1:((p-1)*(p-2)/2));
Nr=[Nrv,Nrs,Nrb];
