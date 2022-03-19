function F=force_integration_vol(etpl,ed,coord,~,~)
% Force integration over the volume of elements.
% Currently only a force integration for the Newmann boundary conditions
% this function needs to be expanded to included all force types
%
% Input(s):
% etpl       - Element tolopogy struture (see seed_mesh.m)
% ed         - Degrees of freedom steering matrix
% coord      - Element coordinates (see seed_mesh.m)
%
% Ouput(s):
% F - Force vector from Newmann/Dirichlet boundary condition

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

p_max=max(etpl.poly(etpl.tree(:,1)==1,2)); nD = 2;                         % Maximum polynomial order
[ngp,vwp,dNr,Nr]=dNr_vol_gauss(p_max);                                     % Shape functions for all elements
F=zeros(max(ed(:)),1);
nels = size(etpl.mat,1);                                                   % Number of elements in the entire mesh
for nel = 1:nels                                                           % Integral loop over all elements
    if etpl.tree(nel,1)==1
        loc_p=etpl.poly(nel,2);                                            % Local polynomial order
        nov=nov_calc(loc_p)*2;                                             % Number of variables
        f=zeros(nov,1);                                                    % Local stiffness matrix
        JT=dNr(:,1:3)*coord(etpl.mat(nel,:),:);                            % Jocabian
        dNrt=choose_dNr(dNr,loc_p,p_max);                                  % Shape function derivatives for this element
        Nrt=choose_Nr(Nr,loc_p,p_max);                                     % Shape functions for this element
        B=zeros(3,nov);                                                    % Assigning space for derivative matrix
        N=zeros(2,nov);                                                    % Assigning space for shapefunction matrix
        for gp = 1:ngp                                                     % Gauss point loop
            indx=2*(gp-1)+(1:nD);                                          % Jacobian and shape function derivative index
            dNx=JT(indx,:)\dNrt(indx,:);                                   % Shape function derivative
            Nx=Nrt(gp,:);                                                  % Shape functions
            B([1 3],1:2:end)=dNx;                                          % Derivative matrix
            B([3 2],2:2:end)=dNx;                                          %
            N(1,1:2:end)=Nx;                                               % Shape functions
            N(2,2:2:end)=Nx;                                               %
            detJ=det(JT(indx,:));                                          % Jacobian
            x=N(1,1:2:6)*coord(etpl.mat(nel,:),1);
            y=N(1,1:2:6)*coord(etpl.mat(nel,:),2);
            f=f+N'*BodyForce(x,y)*detJ*vwp(gp);                               % Local stiffness matrix summation
        end
        DOF=ed(nel,ed(nel,:)~=0);
        F(DOF)=F(DOF)+f;
    end
end
