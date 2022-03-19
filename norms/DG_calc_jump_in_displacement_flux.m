function Er=DG_calc_jump_in_displacement_flux(coord,etpl,etpl_face,ed,u,E,v)
% Error for the jump in displacement
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
% Er     - Error in displacement jump for internal face in the mesh

%  Copyright (C) 2017 Robert Bird
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

Er=zeros(size(etpl.mat,1),1);                                              % Error storage vector
etpl_face=etpl_face(etpl_face(:,end)==1,:);                                % All of the internal faces of the mesh
n_faces=size(etpl_face,1);                                                 % Number of faces
el_p=etpl_face(:,1); el_n=etpl_face(:,2);                                  % Positive and negative element numbers
f_p=etpl_face(:,3);                                                        % Positive and negative local face numbers
nx=etpl_face(:,5);   ny=etpl_face(:,6);                                    % Normals away from small element | (face length only for negative elements)
De=d_mat(E,v);
for fn = 1:n_faces                                                         % Surface integral loop
    if el_n(fn)~=0                                                         % If internal face
        if etpl.tree(el_p(fn),1)==1 && etpl.tree(el_n(fn),1)==1
            co_pve=coord(etpl.mat(el_p(fn),:),:);                          % Positive element coordinates
            co_nve=coord(etpl.mat(el_n(fn),:),:);                          % Negative element coordinates
            pos_loc_p=etpl.poly(el_p(fn),2);                               % Positive element polynomial order
            neg_loc_p=etpl.poly(el_n(fn),2);                               % Negative element polynomial order
            nov_n=nov_calc(neg_loc_p)*2;                                   % positive element number of variables
            p_f=f_p(fn);                                                   % Positive face
            [~,h_small] = face_flag(etpl_face(fn,:),co_pve,co_nve);        % Length of positive
            p_max=max([pos_loc_p,neg_loc_p]);
            [xsi_l,eta_l,xsi_s,eta_s,swp,ngp]=...                          % Gauss point generation
            Gauss_point_positions(p_max,el_p(fn),el_n(fn),coord,etpl.mat,p_f); %
            [~,Nr]=dNr_surf_gauss_mixed(p_max,xsi_s,eta_s,xsi_l,eta_l);    % Shape function
            Nnm=zeros(2,nov_n);
            nov_p=nov_calc(pos_loc_p)*2;                                   % positive element number of variables
            Npm=zeros(2,nov_p);
            for gp = 1:ngp                                                 % Gauss point loop
                indx_sNr = gp;                                             % Shape function index
                W = swp(gp)*h_small/2;
                
                Nrp=choose_Nr(Nr.p,pos_loc_p,p_max); Np = Nrp(indx_sNr,:); % Positive element shape functions
                Npm(1,1:2:end) = Np; Npm(2,2:2:end) = Np;                  % Shape function matrix positive face
                up=Npm*u(ed(el_p(fn),ed(el_p(fn),:)~=0));                  % Postive element displacement
                u_p=up'*[nx(fn);ny(fn)]; 	                               % up multiplyed its normal
                
                Nrn=choose_Nr(Nr.n,neg_loc_p,p_max); Nn = Nrn(indx_sNr,:); % Negative element shape functions
                Nnm(1,1:2:end) = Nn; Nnm(2,2:2:end) = Nn;                  % Shape function matrix negative face
                un=Nnm*u(ed(el_n(fn),ed(el_n(fn),:)~=0));                  % Displacement from the negative element
                u_n=-un'*[nx(fn);ny(fn)];
                
                [pen,~]=pen_func(co_pve,co_nve,h_small,pos_loc_p,neg_loc_p,E,v,nx(fn),ny(fn));
                u2 =(u_p+u_n)^2;
                
                Er(el_p(fn))=Er(el_p(fn))+pen*u2*W;       % Positive error summation
                Er(el_n(fn))=Er(el_n(fn))+pen*u2*W;       % Negative error summation
            end
        end
    end
end