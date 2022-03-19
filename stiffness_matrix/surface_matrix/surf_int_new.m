function [k]=surf_int_new(coord,etpl,etpl_face,E,v,ed)
% SIPG surface integral
%
% Input(s):
% coord      - Element coordinates (see seed_mesh.m)
% etpl       - Element tolopogy structure (see seed_mesh.m)
% etpl_face  - Element face tolopogy structure (see seed_mesh.m)
% E          - Youngs modulus
% v          - Poisson's ratio
% ed         - Current steering matrix
%
% Ouput(s):
% k          - Structure containing all the components of the surface integral
%            .v1, .v2, .v3, .v4 Vectors of the surface integral for surface matrices (1,2,3,4)
%            .r1, .r2, .r3, .r4 Row steering vector of k(.v1,.v2,.v3.v4) -> global stiffness matrix
%            .c1, .c2, .c3, .c4 Column steering vector of k(.v1,.v2,.v3.v4) -> global stiffness matrix

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

etpl_face = etpl_face(etpl_face(:,end)==1,:);                              % All interior face information
nD = 2;                                                                    % number of dimensions
n_faces = size(etpl_face,1);                                               % Number of faces to integrate over
tndof_pn_np=tndof_pn_np_f (etpl.poly,nD,etpl_face);                        % Total number of degrees of freedom for surface matrix 2 and 3
tndof_pp=tndof_pp_f       (etpl.poly,nD,etpl_face);                        % Total number of degrees of freedom for surface matrix 1
tndof_nn=tndof_nn_f       (etpl.poly,nD,etpl_face);                        % Total number of degrees of freedom for surface matrix 4
k.v1=zeros(tndof_pp,1);    k.v2=zeros(tndof_pn_np,1);                      % Zeroed vectors for space for local surface stiffness values
k.v3=zeros(tndof_pn_np,1); k.v4=zeros(tndof_nn,1);                         %
k.r1=zeros(tndof_pp,1);   k.c1=zeros(tndof_pp,1);    kloc1=0;              % Zeroed vectors for row and column steering matrix: Stiffness matrix 1
k.r2=zeros(tndof_pn_np,1);k.c2=zeros(tndof_pn_np,1); kloc2=0;              % Zeroed vectors for row and column steering matrix: Stiffness matrix 2
k.r3=zeros(tndof_pn_np,1);k.c3=zeros(tndof_pn_np,1); kloc3=0;              % Zeroed vectors for row and column steering matrix: Stiffness matrix 3
k.r4=zeros(tndof_nn,1);   k.c4=zeros(tndof_nn,1);    kloc4=0;              % Zeroed vectors for row and column steering matrix: Stiffness matrix 4
nx = etpl_face(:,5); ny = etpl_face(:,6);                                  % Outward normal components
el_p = etpl_face(:,1); el_n = etpl_face(:,2);                              % Positive and negative element numbers
f_p  = etpl_face(:,3); f_n  = etpl_face(:,4);                              % Positive element face numbers
De=d_mat(E,v);                                                             % De = (E/((1+v)*(1-(2*v)))).*[1-v v 0;v 1-v 0;0 0 (1-(2*v))/2];                 % Hookian stiffness matrix
for fn = 1:n_faces                                                         % Surface integral loop
    if etpl.tree(el_p(fn),1)==1 && etpl.tree(el_n(fn),1)==1
        co_pve=coord(etpl.mat(el_p(fn),:),:);                              % Positive element coordinates
        co_nve=coord(etpl.mat(el_n(fn),:),:);                              % Negative element coordinates
        n = [nx(fn);ny(fn)];                                               % Outward normal of face
        pos_loc_p=etpl.poly(el_p(fn),2);                                   % positive element polynomial
        neg_loc_p=etpl.poly(el_n(fn),2);                                   % negative element polynomial
        nov_p=nov_calc(pos_loc_p)*2;                                       % positive element number of variables
        nov_n=nov_calc(neg_loc_p)*2;                                       % positive element number of variables
        
        ke_1 = zeros(nov_p);                                               % Assigning space for temporary local surface matrices: matrix 1
        ke_2 = zeros(nov_p,nov_n);                                         % matrix 2
        ke_3 = zeros(nov_n,nov_p);                                         % matrix 3
        ke_4 = zeros(nov_n);                                               % matrix 4
        p_f=f_p(fn);                                                       % Positive local face number
        [~,h_small] = face_flag(etpl_face(fn,:),co_pve,co_nve);            % Smallest face
        
        [~,h_pos] = face_flag([el_p(fn) 0 f_p(fn) 0],co_pve,co_pve);            % Smallest face
        [~,h_neg] = face_flag([el_n(fn) 0 f_n(fn) 0],co_nve,co_nve);            % Smallest face
        
        if h_neg<h_pos
           fprintf('Error in mesh refinement\n');
        end
        max_p=max([pos_loc_p,neg_loc_p]);
        [xsi_l,eta_l,xsi_s,eta_s,swp,ngp]=Gauss_point_positions(max_p,el_p(fn),el_n(fn),coord,etpl.mat,p_f); % Gauss point positions
        [dNr,Nr]=dNr_surf_gauss_mixed(max_p,xsi_s,eta_s,xsi_l,eta_l);      % Shape functions and their derivatives
        [pen,~]=pen_func(co_pve,co_nve,h_small,pos_loc_p,neg_loc_p,E,v,nx(fn),ny(fn)); % Penalty function
        Bp=zeros(3,nov_p);Npm=zeros(2,nov_p);
        Bn=zeros(3,nov_n);Nnm=zeros(2,nov_n);
        nm = [nx(fn) 0 ny(fn);0 ny(fn) nx(fn)];                            % Normal matrix
        for gp = 1:ngp                                                     % Gauss point loop
            indx_dNr =((gp-1)*nD)+(1:nD);                                  % Shape function derivative index
            indx_sNr = gp;                                                 % Shape function index
            Jp = dNr.p(indx_dNr,1:3)*co_pve;                               % Jacobian for positive element
            Jn = dNr.n(indx_dNr,1:3)*co_nve;                               % Jacobian for negative element
            dNrp=choose_dNr(dNr.p,pos_loc_p,max_p); dNx_p=Jp\dNrp(indx_dNr,:); % Global shape function derivatives positive
            dNrn=choose_dNr(dNr.n,neg_loc_p,max_p); dNx_n=Jn\dNrn(indx_dNr,:); % Global shape function derivatives negative
            Nrp=choose_Nr(Nr.p,pos_loc_p,max_p); Np = Nrp(indx_sNr,:);     % Global shape functions positive
            Nrn=choose_Nr(Nr.n,neg_loc_p,max_p); Nn = Nrn(indx_sNr,:);     % Global shape functions negative
            Bp([1 3],1:2:end)=dNx_p; Bp([3 2],2:2:end)=dNx_p;              % Derivative matrix positive face
            Bn([1 3],1:2:end)=dNx_n; Bn([3 2],2:2:end)=dNx_n;              % Derivative matrix negative face
            Npm(1,1:2:end) = Np; Npm(2,2:2:end) = Np;                      % Shape function matrix positive face
            Nnm(1,1:2:end) = Nn; Nnm(2,2:2:end) = Nn;                      % Shape function matrix negative face
            W = swp(gp)*h_small/2;                                         % Integral weight

            ke_1 = ke_1-((Bp'*De'*nm'*Npm)/2)*W-((Npm'*nm*De*Bp)/2)*W+pen*(Npm'*Npm)*W; % Computing the surface stiffness matrix 1
            ke_2 = ke_2+((Bp'*De'*nm'*Nnm)/2)*W-((Npm'*nm*De*Bn)/2)*W-pen*(Npm'*Nnm)*W; % Computing the surface stiffness matrix 2
            ke_3 = ke_3-((Bn'*De'*nm'*Npm)/2)*W+((Nnm'*nm*De*Bp)/2)*W-pen*(Nnm'*Npm)*W; % Computing the surface stiffness matrix 3
            ke_4 = ke_4+((Bn'*De'*nm'*Nnm)/2)*W+((Nnm'*nm*De*Bn)/2)*W+pen*(Nnm'*Nnm)*W; % Computing the surface stiffness matrix 4
        end
        neDoF_1 = ((nov_p)^2);                                             % Number of terms in local stiffness matrix 1
        neDoF_2 = ((nov_p*nov_n));                                         % Number of terms in local stiffness matrix 2
        neDoF_3 = ((nov_p*nov_n));                                         % Number of terms in local stiffness matrix 3
        neDoF_4 = ((nov_n)^2);                                             % Number of terms in local stiffness matrix 4
        kloc1 = max(kloc1)+(1:neDoF_1);                                    % Steering matrix 1 index
        kloc2 = max(kloc2)+(1:neDoF_2);                                    % Steering matrix 2 index
        kloc3 = max(kloc3)+(1:neDoF_3);                                    % Steering matrix 3 index
        kloc4 = max(kloc4)+(1:neDoF_4);                                    % Steering matrix 4 index
        k.v1(kloc1) = reshape(ke_1,neDoF_1,1);                             % Storage of local stiffness matrix 1
        k.v2(kloc2) = reshape(ke_2,neDoF_2,1);                             % Storage of local stiffness matrix 2
        k.v3(kloc3) = reshape(ke_3,neDoF_3,1);                             % Storage of local stiffness matrix 3
        k.v4(kloc4) = reshape(ke_4,neDoF_4,1);                             % Storage of local stiffness matrix 4
        ed_r=ed(el_p(fn),1:nov_p);                                         % Steering vector rows for matrix 1
        k.r1(kloc1)=reshape(ed_r.'*ones(1,nov_p),neDoF_1,1);               %
        ed_c=ed(el_p(fn),1:nov_p);                                         % Steering vector columns for matrix 1
        k.c1(kloc1)=reshape(ones(nov_p,1)*ed_c  ,neDoF_1,1);               %
        ed_r=ed(el_p(fn),1:nov_p);                                         % Steering vector rows for matrix 2
        k.r2(kloc2)=reshape(ed_r.'*ones(1,nov_n),neDoF_2,1);               %
        ed_c=ed(el_n(fn),1:nov_n);                                         % Steering vector columns for matrix 2
        k.c2(kloc2)=reshape(ones(nov_p,1)*ed_c  ,neDoF_2,1);               %
        ed_r=ed(el_n(fn),1:nov_n);                                         % Steering vector rows for matrix 3
        k.r3(kloc3)=reshape(ed_r.'*ones(1,nov_p),neDoF_3,1);               %
        ed_c=ed(el_p(fn),1:nov_p);                                         % Steering vector columns for matrix 3
        k.c3(kloc3)=reshape(ones(nov_n,1)*ed_c  ,neDoF_3,1);               %
        ed_r=ed(el_n(fn),1:nov_n);                                         % Steering vector rows for matrix 4
        k.r4(kloc4)=reshape(ed_r.'*ones(1,nov_n),neDoF_4,1);               %
        ed_c=ed(el_n(fn),1:nov_n);                                         % Steering vector columns for matrix 4
        k.c4(kloc4)=reshape(ones(nov_n,1)*ed_c  ,neDoF_4,1);               %
    end
end