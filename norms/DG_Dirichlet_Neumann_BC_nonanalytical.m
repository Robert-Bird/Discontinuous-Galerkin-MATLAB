function Er=DG_Dirichlet_Neumann_BC_nonanalytical(coord,etpl,etpl_face,ed,u,E,v)
% Estimated error for the Dirichlet / Neumann boundary condition
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
% Er     - Error for all faces on the Dirichlet / Neumann BC

%  Copyright (C) 2017 Robert Bird
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

Er=zeros(size(etpl.mat,1),1);                                              % Error storage vector
p_max=max(etpl.poly(:,2)); nD = 2;                                         % Maximum polynomial order for the external faces
[swp,~,Nr,ngp]=dNr_surf_gauss(p_max);                                      % Shape functions, derivatives and guass points
etpl_face=etpl_face(abs(etpl_face(:,end))==4,:);                           % All Dirichlet BCs
n_faces=size(etpl_face,1);                                                 % Number of faces
el_p=etpl_face(:,1);                                                       % Positive and negative element numbers
f_p=etpl_face(:,3);                                                        % Positive and negative local face numbers
nx=etpl_face(:,5);   ny=etpl_face(:,6);                                    % Normals away from small element | (face length only for negative elements)
for fn = 1:n_faces                                                         % Surface integral loop
    if etpl.tree(el_p(fn),1)==1
        co_pve=coord(etpl.mat(el_p(fn),:),:);
        p_f=f_p(fn);                                                       % External element local face number
        pos_loc_p=etpl.poly(el_p(fn),2);                                   % External element polynomial order
        nov_p=nov_calc(pos_loc_p)*2;                                       % positive element number of variables
        Npm=zeros(2,nov_p);
        h_small=face_size_cal(p_f,coord(etpl.mat(el_p(fn),:),:));
        for gp = 1:ngp                                                     % Gauss point loop
            indx_sNr = gp+((p_f-1)*ngp);                                   % Shape function index
            Nrp=choose_Nr(Nr,pos_loc_p,p_max);   Np = Nrp(indx_sNr,:);     % Positive element shape functions
            Npm(1,1:2:end) = Np; Npm(2,2:2:end) = Np;                      % Shape function matrix positive face
            W = swp(gp)*h_small/2;                                         % Integral weight
            u_h=Npm*u(ed(el_p(fn),ed(el_p(fn),:)~=0));                     % Postive element displacement
            
            u2=(u_h'*[nx(fn);ny(fn)])^2;
            
            
            [pen,~]=pen_func(co_pve,co_pve,h_small,pos_loc_p,pos_loc_p,E,v,nx(fn),ny(fn));
            Er(el_p(fn))=Er(el_p(fn))+pen*u2*W;
        end
    end
end

