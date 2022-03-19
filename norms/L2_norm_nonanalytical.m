function [L2] = L2_norm_nonanalytical(coord,etpl,ed,u)
% L2 norm calculation for problem without analytical solution.
% The reference solution is the zero solution.
%
% Input(s):
% coord      - Coordinates of all nodes in the mesh
% etpl       - Element topology structure
% ed         - Degree of freedom steering matrix
% u          - Displacment solution
%
% Ouput(s):
% L2 - L2 norm

%  Copyright (C) 2017 Robert Bird
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

p_max=max(etpl.poly(etpl.tree(:,1)==1,2)); nD = 2;                         % Maximum polynomial order
[ngp,vwp,dNr,Nr]=dNr_vol_gauss(p_max);                                     % Shape functions for all elements
L2  = 0;                                                                   % L2 norm
nels = size(etpl.mat,1);                                                   % number of elements in the mesh (active and deactivated)
for nel = 1:nels                                                           % Integral loop over all elements
    if etpl.tree(nel,1)==1                                                 % Check to make sure the element is active
        loc_p=etpl.poly(nel,2);                                            % Local polynomial order
        nov=nov_calc(loc_p)*2;                                             % Number of variables
        JT=dNr(:,1:3)*coord(etpl.mat(nel,:),:);                            % Jacobian
        Nrt=choose_Nr(Nr,loc_p,p_max);                                     % Shape functions for this element
        N=zeros(2,nov);                                                    % Assigning space for shapefunction matrix
        for gp = 1:ngp                                                     % Gauss point loop
            indx=2*(gp-1)+(1:nD);                                          % Jacobian index
            Nx=Nrt(gp,:);                                                  % Shape functions
            N(1,1:2:end)=Nx;                                               % Shape function matrix
            N(2,2:2:end)=Nx;                                               %
            X=Nx(1:3)*coord(etpl.mat(nel,:),1);                            % X-coordinate for analytical expression
            Y=Nx(1:3)*coord(etpl.mat(nel,:),2);                            % Y-coordinate for analytical expression
            U=[0 0];
            u_analytical(1)=U(1);                                          % Analytical expression for u and v e.g. sin(pi*x)*sin(pi*y)
            u_analytical(2)=U(2);                                          %
            u_computed=N*u(ed(nel,ed(nel,:)~=0));                          % Computed solution
            norm_squared=(u_computed(1)-u_analytical(1))^2+(u_computed(2)-u_analytical(2))^2; % the norm square at a point
            detJ=det(JT(indx,:));                                          % Jacobian determinant
            L2=L2+norm_squared*detJ*vwp(gp);                               % L2 norm calculation
        end
    end
end
L2=sqrt(L2);                                                               % Square rooting the calculation to give the completed norm
