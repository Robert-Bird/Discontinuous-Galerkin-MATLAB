function Er=Error_calc_Neumann_BC_analytical(coord,etpl,etpl_face,ed,u,E,v)
% Estimated error for the Neumann boundary condition
%
% Input(s):
% coord      - Coordinates of all nodes in the mesh
% etpl       - Element topology structure
% etpl_face  - Element face topology matrix
% ed         - Degree of freedom steering matrix
% u          - Displacment solution
% E          - Young's modulus
% v          - Poisson's ratio
%
% Ouput(s):
% Er     - Error for all faces on the Neumann BC

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

Er=zeros(size(etpl.mat,1),1);    nD=2;                                     % Error storage vector | and number of dimensions
p_max=max(etpl.poly(:,2));                                                 % Maximum polynomial order for the external faces
[swp,dNr,Nr,ngp]=dNr_surf_gauss(p_max);                                     % Shape functions, derivatives and guass points
etpl_face=etpl_face(abs(etpl_face(:,end))==3,:);                           % Internal etpl_face
n_faces=size(etpl_face,1);                                                 % Number of faces
el_p=etpl_face(:,1);                                                       % Positive and negative element numbers
f_p=etpl_face(:,3);                                                        % Positive and negative local face numbers
nx=etpl_face(:,5);   ny=etpl_face(:,6);                                    % Normals away from small element
De=d_mat(E,v);                                                             % Hookian stiffness matrix
for fn = 1:n_faces                                                         % Surface integral loop
    nm = [nx(fn) 0 ny(fn);0 ny(fn) nx(fn)];                                % Normal matrix
    co_pve=coord(etpl.mat(el_p(fn),:),:);                                  % Coordinates of negative element
    pos_loc_p=etpl.poly(el_p(fn),2);                                       % positive element polynomial
    nov_p=nov_calc(pos_loc_p)*2;                                           % positive element number of variables
    Bp=zeros(3,nov_p);
    p_f=f_p(fn);                                                           % Positive face number
    h_small=face_size_cal(p_f,coord(etpl.mat(el_p(fn),:),:));              % Finding the smallest face
    for gp = 1:ngp                                                         % Surface gauss point loop
        indx_dNr = ((gp-1)*nD)+((p_f-1)*nD*ngp)+(1:nD);                    % Positive shape function derivative index
        indx_sNr = gp+((p_f-1)*ngp);                                       % Shape function index
        Nrp=choose_Nr(Nr,pos_loc_p,p_max); Np = Nrp(indx_sNr,:);           % Positive element shape functions
        Jp = dNr(indx_dNr,1:3)*co_pve;                                     % Jacobian of positive element
        dNrp=choose_dNr(dNr,pos_loc_p,p_max); dNx_p=Jp\dNrp(indx_dNr,:);   % Positive global shape function derivative
        Bp([1 3],1:2:end)=dNx_p; Bp([3 2],2:2:end)=dNx_p;                  % Derivative matrix positive face
        W = swp(gp)*h_small/2;                                             % Integral weight
        up=u(ed(el_p(fn),ed(el_p(fn),:)~=0));                              % Positive element values for u
        s_p_norm=nm*De*Bp*up;                                              % Finite element approximation
        
        
        if etpl_face(fn,end)==3
            u_f=(s_p_norm(1))^2 + (s_p_norm(2))^2;
        elseif etpl_face(fn,end)~=3
            X=(Np(1:3)*coord(etpl.mat(el_p(fn),:),1));                     % X coordinate
            Y=(Np(1:3)*coord(etpl.mat(el_p(fn),:),2));                     % Y coordinate
            s_n_norm=[nx(fn) 0 ny(fn);0 ny(fn) nx(fn)]*BodyStress(X,Y);        % Applied Neumann BC
            u_f=(s_p_norm(1)-s_n_norm(1))^2 + (s_p_norm(2)-s_n_norm(2))^2; % Difference in Neuman BC
        end
        Er(el_p(fn))=Er(el_p(fn))+(h_small/(pos_loc_p))*u_f*W;             % Error to positive element
    end 
end

end

