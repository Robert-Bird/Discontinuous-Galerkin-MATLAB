function Er=Error_calc_jump_in_stress_flux(coord,etpl,etpl_face,ed,u,E,v)
% Error for the jump in stress in internal faces
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
% Er     - Error in stress jump for every internal face in the mesh

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

Er=zeros(size(etpl.mat,1),1);        nD=2;                                 % Error storage vector | and number of dimensions
etpl_face=etpl_face(etpl_face(:,end)==1,:);                                % Internal etpl_face
n_faces=size(etpl_face,1);                                                 % Number of faces
el_p=etpl_face(:,1); el_n=etpl_face(:,2);                                  % Positive and negative element numbers
f_p=etpl_face(:,3);                                                        % Positive local face numbers
nx=etpl_face(:,5);   ny=etpl_face(:,6);                                    % Normals away from small element
De=d_mat(E,v);                                                             % Hookian stiffness matrix
for fn = 1:n_faces                                                         % Surface integral loop
    if etpl.tree(el_p(fn),1)==1 && etpl.tree(el_n(fn),1)==1
        nm = [nx(fn) 0 ny(fn);0 ny(fn) nx(fn)];                            % Normal matrix
        co_pve=coord(etpl.mat(el_p(fn),:),:);                              % Coordinates of negative element
        co_nve=coord(etpl.mat(el_n(fn),:),:);                              % Coordinates of negative element
        pos_loc_p=etpl.poly(el_p(fn),2);                                   % positive element polynomial
        neg_loc_p=etpl.poly(el_n(fn),2);                                   % negative element polynomial
        nov_p=nov_calc(pos_loc_p)*2;                                       % positive element number of variables
        nov_n=nov_calc(neg_loc_p)*2;                                       % positive element number of variables
        Bp=zeros(3,nov_p);
        Bn=zeros(3,nov_n);
        p_f=f_p(fn);                                                       % Positive face number
        p_max=max([pos_loc_p,neg_loc_p]);
        [~,h_small] = face_flag(etpl_face(fn,:),co_pve,co_nve);            % Finding the smallest face
        [xsi_l,eta_l,xsi_s,eta_s,swp,ngp]=Gauss_point_positions(p_max,el_p(fn),el_n(fn),coord,etpl.mat,p_f); % Small and large face Gauss point generation
        [dNr,~]=dNr_surf_gauss_mixed(p_max,xsi_s,eta_s,xsi_l,eta_l);       % Local shape function derivatives for small and large
        for gp = 1:ngp                                                     % Surface gauss point loop
            indx_dNr =((gp-1)*nD)+(1:nD);                                  % Shape function index
            Jp = dNr.p(indx_dNr,1:3)*co_pve;                               % Jacobian of positive element
            Jn = dNr.n(indx_dNr,1:3)*co_nve;                               % Jacobian of negative element
            dNrp=choose_dNr(dNr.p,pos_loc_p,p_max); dNx_p=Jp\dNrp(indx_dNr,:); % Positive global shape function derivative
            dNrn=choose_dNr(dNr.n,neg_loc_p,p_max); dNx_n=Jn\dNrn(indx_dNr,:); % Negative global shape function derivative
            Bp([1 3],1:2:end)=dNx_p; Bp([3 2],2:2:end)=dNx_p;              % Derivative matrix positive face
            Bn([1 3],1:2:end)=dNx_n; Bn([3 2],2:2:end)=dNx_n;              % Derivative matrix negative face
            max_p=max([pos_loc_p,neg_loc_p]);                              % Maximum polynomial order
            W = swp(gp)*h_small/2;                                         % Integral weight
            up=u(ed(el_p(fn),ed(el_p(fn),:)~=0));                          % Positive element values for u
            un=u(ed(el_n(fn),ed(el_n(fn),:)~=0));                          % Negative element values for u
            s_p_norm=nm*De*Bp*up;
            s_n_norm=nm*De*Bn*un;
            u_f=(s_p_norm(1)-s_n_norm(1))^2 + (s_p_norm(2)-s_n_norm(2))^2;
            Er(el_p(fn))=Er(el_p(fn))+0.5*(h_small/(max_p))*u_f*W;         % Error to positive element
            Er(el_n(fn))=Er(el_n(fn))+0.5*(h_small/(max_p))*u_f*W;         % Error to negative element
        end
    end
end
