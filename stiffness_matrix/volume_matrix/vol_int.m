function [kval,krow,kcol,ed_new] = vol_int(coord,etpl,E,v)
% SIPG volumetric integral
%
% Input(s):
% coord      - Element coordinates (see seed_mesh.m)
% etpl       - Element tolopogy struture (see seed_mesh.m)
% E          - Youngs modulus
% v          - Poisson's ratio
%
% Ouput(s):
% kval             - Vector of all local element volumetric stiffness components
% krow             - Row steering vector of kval -> global stiffness matrix
% kcol             - Column steering vector of kval -> global stiffness matrix
% ed_new           - Updated steering matrix components

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

p_max=max(etpl.poly(etpl.tree(:,1)==1,2)); nD = 2;                         % Maximum polynomial order
[ngp,vwp,dNr,Nr]=dNr_vol_gauss(p_max);                                     % Shape functions for all elements
ed_new=zeros(size(etpl.mat,1),(nov_calc(p_max))*nD);                       % Zeroes matrix for degrees of freedom
tndof=tndof_sum(etpl.poly(etpl.tree(:,1)==1,:),2);                         % Total number of degrees of freedom
krow  = zeros(tndof,1); kcol=krow; kval=krow; kloc=0;                      % Zeroed steering vectors
De=d_mat(E,v);                                                             % Hookian stiffness matrix
nels = size(etpl.mat,1);                                                   % Number of elements in the entire mesh
for nel = 1:nels                                                           % Integral loop over all elements
    if etpl.tree(nel,1)==1                                                                                             
        loc_p=etpl.poly(nel,2);                                            % Local polynomial order
        nov=nov_calc(loc_p)*nD;                                            % Number of variables
        ke=zeros(nov);                                                     % Local stiffness matrix
        JT=dNr(:,1:3)*coord(etpl.mat(nel,:),:);                            % Jocabian
        dNrt=choose_dNr(dNr,loc_p,p_max);                                  % LOCAL Shape function derivatives for this element
        Nrt=choose_Nr(Nr,loc_p,p_max);                                     % Shape functions for this element
        B=zeros(3,nov);                                                    % Assigning space for derivative matrix
        N=zeros(nD,nov);                                                   % Assigning space for shapefunction matrix
        for gp = 1:ngp                                                     % Gauss point loop
            indx=nD*(gp-1)+(1:nD);                                         % Jacobian and shape function derivative index
            dNx=JT(indx,:)\dNrt(indx,:);                                   % Shape function derivative
            Nx=Nrt(gp,:);                                                  % Shape functions
            B([1 3],1:2:end)=dNx;                                          % Derivative matrix
            B([3 2],2:2:end)=dNx;                                          %
            N(1,1:2:end)=Nx;                                               % Shape functions
            N(2,2:2:end)=Nx;                                               %
            detJ=det(JT(indx,:));                                          % Jacobian
            if detJ<=0
                fprintf('Volume <0\n');                                    % Checking if the element is inverted
            end
            ke=ke+B.'*De*B*detJ*vwp(gp);                                   % Local stiffness matrix summation
        end
        ed_new(nel,1:nov)=(1:nov)+max(ed_new(:));                          % Creating degrees of freedom for current element
        ed_t = ed_new(nel,1:nov);                                          % degrees of freedom for current element
        kloc=(1:nov^2)+max(kloc);                                          % Steering vector index
        neDoF=(nov)^2;                                                     % Number of terms in local matrix
        krow(kloc)=reshape(ed_t.'*ones(1,nov),neDoF,1);                    % Steering vector row values
        kcol(kloc)=reshape(ones(nov,1)*ed_t  ,neDoF,1);                    % Steering vector column values
        kval(kloc)=reshape(ke,neDoF,1);                                    % Vector of local stiffness values
    end
end
