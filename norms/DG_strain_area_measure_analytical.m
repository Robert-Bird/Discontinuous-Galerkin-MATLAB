function Er=DG_strain_area_measure_analytical(coord,etpl,ed,u,E,v)
% Estimated error for the approximation of the strong form
%
% Input(s):
% coord      - Coordinates of all nodes in the mesh
% etpl       - Element topology structure
% ed         - Degree of freedom steering matrix
% u          - Displacment solution
% E          - Young's modulus
% v          - Poisson's ratio
%
% Ouput(s):
% Er     - Strong form error for every element within the mesh, squared

%  Copyright (C) 2017 Robert Bird
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

Er=zeros(size(etpl.mat,1),1);
p_max=max(etpl.poly(etpl.tree(:,1)==1,2)); nD = 2;                         % Maximum polynomial order
[ngp,vwp,dNr,Nr]=dNr_vol_gauss(p_max);                                     % Shape functions for all elements
nels = size(etpl.mat,1);                                                   % Number of elements in the entire mesh
C=inv(d_mat(E,v));
for nel = 1:nels                                                           % Integral loop over all elements
    if etpl.tree(nel,1)==1
        loc_p=etpl.poly(nel,2);                                            % Local polynomial order
        nov=nov_calc(loc_p)*2;                                             % Number of variables
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
            
            Strain_h=B*u(ed(nel,ed(nel,:)~=0));
            
            x=N(1,1:2:6)*coord(etpl.mat(nel,:),1);
            y=N(1,1:2:6)*coord(etpl.mat(nel,:),2);
            strain_analytical=C*BodyStress(x,y); 
            strain_diff=norm([Strain_h-strain_analytical])^2;
            Er(nel)=Er(nel)+strain_diff*detJ*vwp(gp);                      % Local stiffness matrix summation
        end
    end
end

