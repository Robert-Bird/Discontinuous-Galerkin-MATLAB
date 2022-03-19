function [k]=surf_dirichlet_newmann(coord,etpl,etpl_face,E,v,ed)
% SIPG homogeneous Dirichlet / Neumann boundary condition
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
%            .v1 Vectors of the surface integral for the Dirichlet / Neumann BC
%            .r1 Row steering vector of kv.1 -> global stiffness matrix
%            .c1 Column steering vector of k.v1 -> global stiffness matrix

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

etpl_face = etpl_face(abs(etpl_face(:,end))==4,:);                         % Face information for external elements
p_max=max(etpl.poly(:,2)); nD = 2;                                         % Maximum polynomial order for the external faces
[swp,dNr,Nr,ngp]=dNr_surf_gauss(p_max);                                    % Shape functions, derivatives and guass points
n_faces = size(etpl_face,1);                                               % Number of external faces
tndof_pp=tndof_pp_f(etpl.poly,nD,etpl_face);                               % Total number of degrees of freedom for all local stiffness matrices
k.v1 = zeros(tndof_pp,1);                                                  % Assigning space for local surface integrals
k.r1 = zeros(tndof_pp,1);k.c1 = zeros(tndof_pp,1); kloc1=0;                % Assigning space for steering vectors | index counter
nx   = etpl_face(:,5);     ny = etpl_face(:,6);                            % Outward normal componenets and lenght of element face.
el_p = etpl_face(:,1);   f_p  = etpl_face(:,3);                            % Positive element number and faces.
De=d_mat(E,v);                                                             % Hookian stiffness matrix
for fn = 1:n_faces                                                         % Integral loop
    if etpl.tree(el_p(fn),1)==1                                            % Does the element exist in the tree
        pos_loc_p=etpl.poly(el_p(fn),2);                                   % Positive element polynomial order
        nov_p=nov_calc(pos_loc_p)*2;                                       % Positive element number of variables
        ke_1 = zeros(nov_p);                                               % Zeroing temporary space for local surface integrals
        p_f=f_p(fn);                                                       % Positive element face number
        Bp=zeros(3,nov_p);Npm=zeros(2,nov_p);
        nm = [nx(fn) 0 ny(fn);0 ny(fn) nx(fn)];                            % Outward normal matrix
        h=face_size_cal(p_f,coord(etpl.mat(el_p(fn),:),:));                % Length of the face
        for gp = 1:ngp
            indx_dNr = ((gp-1)*nD)+((p_f-1)*nD*ngp)+(1:nD);                % Positive shape function derivative index
            indx_sNr = gp+((p_f-1)*ngp);                                   % Positive shape function index
            co_pve   = coord(etpl.mat(el_p(fn),:),:);                      % Coordinates of positive element
            Jp       = dNr(indx_dNr,1:3)*co_pve;                           % Jacobian of positive element
            [pen,~]  = pen_func(co_pve,co_pve,h,pos_loc_p,pos_loc_p,E,v,nx(fn),ny(fn));	% Penalty term
            dNrp     = choose_dNr(dNr,pos_loc_p,p_max); dNx_p=Jp\dNrp(indx_dNr,:);	% Global shape function derivatives
            Nrp      = choose_Nr(Nr,pos_loc_p,p_max); Np = Nrp(indx_sNr,:);	% Local shapen functions
            Bp([1 3],1:2:end) = dNx_p; Bp([3 2],2:2:end)=dNx_p;            % Shape function derivative matrix
            Npm(1,1:2:end)    = Np; Npm(2,2:2:end) = Np;                   % Shape function matrix positive face
            W                 = swp(gp)*h/2;                               % Integral weight
            n2=[nx(fn);ny(fn)];
            K_part=Npm'*n2*(n2')*nm*De*Bp;%
            ke_1 = ke_1 - K_part*W  - K_part'*W + pen*(Npm'*n2*(n2')*Npm)*W;%
        end
        neDoF_1 = ((nov_p)^2);                                             % Number of degrees of freedom in local matrix
        kloc1 = max(kloc1)+(1:neDoF_1);                                    % Steering vector index
        k.v1(kloc1) = reshape(ke_1,neDoF_1,1);                             % Store the newly calculated term
        ed_r=ed(el_p(fn),1:nov_p);                                         % Steering vector rows for matrix
        k.r1(kloc1)=reshape(ed_r.'*ones(1,nov_p),neDoF_1,1);               %
        ed_c=ed(el_p(fn),1:nov_p);                                         % Steering vector columns for matrix
        k.c1(kloc1)=reshape(ones(nov_p,1)*ed_c  ,neDoF_1,1);               %
    end
end