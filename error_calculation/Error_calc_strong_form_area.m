function Er=Error_calc_strong_form_area(coord,etpl,ed,u,E,v)
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

Er=zeros(size(etpl.mat,1),1);                                              % Error storage vector
nels = size(etpl.mat,1);nD=2;                                              % Number of elements | Number of dimensions
p_max=max(etpl.poly(etpl.tree(:,1)==1,2));                                 % Maximum polynomial order
[ngp,vwp,dNr,Nr] = dNr_vol_gauss(p_max);                                   % Shape functions and derivatives
[~,~,ddNr]       = ddNr_vol_gauss(p_max);                                  % Shape functions and 2nd order derivatives xx yy
[~,~,ddNr_xy]    = ddNr_vol_gauss_xy(p_max);                               % Shape functions and 2nd order derivatives xy
De = d_mat(E,v);                                                           % Hookian stiffness matrix
d_11  = De(1,1); d_12  = De(1,2); d_13  = De(1,3);                         % Hookian stiffness matrix components broken down for ease of use
d_21  = De(2,1); d_22  = De(2,2); d_23  = De(2,3);                         %
d_31  = De(3,1); d_32  = De(3,2); d_33  = De(3,3);                         %


for nel = 1:nels                                                           % Element loop
    if etpl.tree(nel,1)==1
        JT=dNr(:,1:3)*coord(etpl.mat(nel,:),:);                            % Jacobian
        loc_p=etpl.poly(nel,2);                                            % Polynomial order of element
        ddNrt=choose_dNr(ddNr,loc_p,p_max);                                % Second order shape functions xx yy
        ddNr_xyt=choose_ddNr_xy(ddNr_xy,loc_p,p_max);                      % Second order shape functions xy
        n1=coord(etpl.mat(nel,1),:);                                       % Node 1 coordinates
        n2=coord(etpl.mat(nel,2),:);                                       % Node 2 coordinates
        n3=coord(etpl.mat(nel,3),:);                                       % Node 3 coordinates
        a=sqrt(((n1(1)-n2(1))^2)+((n1(2)-n2(2))^2));                       % Finding the average element side length
        b=sqrt(((n3(1)-n2(1))^2)+((n3(2)-n2(2))^2));                       %
        c=sqrt(((n1(1)-n3(1))^2)+((n1(2)-n3(2))^2));                       %
        h=2*(a*b*c)/sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c));                 % Inner circle radius of each element (I think, could be maximum)
        for gp = 1:ngp                                                     % Gauss point loop
            indx=2*(gp-1)+(1:nD);                                          % Index for the Jacobian and the second order differential in xx and yy
            detJ=det(JT(indx,:));                                          % Jacobian determinant
            ddNr_all=[ddNrt(indx(1),:);ddNr_xyt(gp,:);ddNrt(indx(2),:)];   % 2nd order derivatives in one matrix
            Jt=inv(JT(indx,:));                                            % Defining the 2nd order Jacobian
            J2=[Jt(1,1)^2,      2*Jt(1,1)*Jt(1,2),                    Jt(1,2)^2;% 2nd order Jacobian matrix
                Jt(1,1)*Jt(2,1),Jt(1,2)*Jt(2,1)+Jt(1,1)*Jt(2,2),Jt(1,2)*Jt(2,2);%
                Jt(2,1)^2,      2*Jt(2,1)*Jt(2,2),                   Jt(2,2)^2];%
            ddNxx=J2*ddNr_all;                                             % Second order derivative components
            ed_t=ed(nel,ed(nel,:)~=0);                                     % Steering matrix for current element
            ddNu_xx=ddNxx(1,:)*u(ed_t(1:2:end));                           % u differentiated wrt xx
            ddNu_xy=ddNxx(2,:)*u(ed_t(1:2:end));                           % u differentiated wrt xy
            ddNu_yy=ddNxx(3,:)*u(ed_t(1:2:end));                           % u differentiated wrt yy
            ddNv_xx=ddNxx(1,:)*u(ed_t(2:2:end));                           % v differentiated wrt xx
            ddNv_xy=ddNxx(2,:)*u(ed_t(2:2:end));                           % v differentiated wrt xy
            ddNv_yy=ddNxx(3,:)*u(ed_t(2:2:end));                           % v differentiated wrt yy
            x=Nr(gp,1:3)*coord(etpl.mat(nel,:),1);
            y=Nr(gp,1:3)*coord(etpl.mat(nel,:),2);
            F=BodyForce(x,y);
            A = (ddNu_xx)*d_11  + (ddNu_xy)*d_31 + (ddNu_xy)*d_13 + (ddNu_yy)*d_33 + ...
                (ddNv_xx)*d_13  + (ddNv_xy)*d_33 + (ddNv_xy)*d_12 + (ddNv_yy)*d_32;
            B = (ddNu_xx)*d_31  + (ddNu_xy)*d_21 + (ddNu_xy)*d_33 + (ddNu_yy)*d_23 + ...
                (ddNv_xx)*d_33  + (ddNv_xy)*d_23 + (ddNv_xy)*d_32 + (ddNv_yy)*d_22;
            
            DIFF= (B+F(2))^2+(A+F(1))^2;%
            Er(nel)=Er(nel)+((h^2)/(loc_p^2))*(DIFF)*detJ*vwp(gp);          % Final summation
        end
    end
end