function F=force_integration_Dirichlet_Neumann(etpl,etpl_face,ed,coord,E,v)
% Dirichlet / Neumann boundary conditions
%
% Input(s):
% etpl       - Element tolopogy struture (see seed_mesh.m)
% etpl_face  - Element face tolopogy struture (see seed_mesh.m)
% ed         - Degrees of freedom steering matrix
% pr         - Pressure
% coord      - Element coordinates (see seed_mesh.m)
% Dirichlet / Neumann  - Dirichlet / Neumann BC function handle
% E          - Young's Modulus
% v          - Poisson's ratio
%
% Ouput(s):
% F - Force vector from Dirichlet / Neumann boundary condition

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

F=zeros(max(max(ed)),1);                                                   % Assigning variable space for the force vector
etpl_face = etpl_face(etpl_face(:,end)==-4,:);                              % Face information for external elements
p_max=max(etpl.poly(:,2)); nD = 2;                                         % Maximum polynomial order for the external faces
[swp,dNr,Nr,ngp]=dNr_surf_gauss(p_max);                                    % Shape functions, derivatives and guass points
n_faces = size(etpl_face,1);                                               % Number of external faces
nx   = etpl_face(:,5);     ny = etpl_face(:,6);                            % Outward normal componenets and lenght of element face.
en = etpl_face(:,1);   f_p  = etpl_face(:,3);                              % Positive element number and faces.
De=d_mat(E,v);                                                             % Hookian stiffness matrix
for i = 1:n_faces                                                          % Integral loop
    pos_loc_p=etpl.poly(en(i),2);                                          % Positive element polynomial order
    nov_p=nov_calc(pos_loc_p)*2;                                           % Positive element number of variables
    Fl = zeros(nov_p,1);                                                   % Zeroing temporary space for local surface integrals
    p_f=f_p(i);                                                            % Positive element face number
    Bp=zeros(3,nov_p);Npm=zeros(2,nov_p);
    nm = [nx(i) 0 ny(i);0 ny(i) nx(i)];                                    % Outward normal matrix
    h=face_size_cal(p_f,coord(etpl.mat(en(i),:),:));                       % Length of the face
    for gp = 1:ngp
        indx_dNr = ((gp-1)*nD)+((p_f-1)*nD*ngp)+(1:nD);                    % Positive shape function derivative index
        indx_sNr = gp+((p_f-1)*ngp);                                       % Positive shape function index
        co_pve   = coord(etpl.mat(en(i),:),:);                             % Coordinates of positive element
        Jp       = dNr(indx_dNr,1:3)*co_pve;                               % Jacobian of positive element
        [pen,~]  = pen_func(co_pve,co_pve,h,pos_loc_p,pos_loc_p,E,v,nx(i),ny(i));	% Penalty term
        dNrp     = choose_dNr(dNr,pos_loc_p,p_max); dNx_p=Jp\dNrp(indx_dNr,:);	% Global shape function derivatives
        Nrp      = choose_Nr(Nr,pos_loc_p,p_max); Np = Nrp(indx_sNr,:);	   % Local shapen functions
        Bp([1 3],1:2:end) = dNx_p; Bp([3 2],2:2:end)=dNx_p;                % Shape function derivative matrix
        Npm(1,1:2:end)    = Np; Npm(2,2:2:end) = Np;                       % Shape function matrix positive face
        W                 = swp(gp)*h/2;                                   % Integral weight
        n2=[nx(i);ny(i)];
        X=(Npm(1,1:2:6)*coord(etpl.mat(en(i),:),1));                       % Global x-position
        Y=(Npm(2,2:2:6)*coord(etpl.mat(en(i),:),2));                       % Global y-position
        u=ImposedDisplacements(X,Y);                                       % Displacement vector at global position (x,y)
        K_part=n2*(n2')*nm*De*Bp;
        Fl = Fl - K_part'*W*u + pen*(Npm'*n2*(n2'))*W*u;%
    end
    ed_s=ed(en(i),1:sum(ed(en(i),:)~=0));                                  % Local degrees of freedom
    F(ed_s)=F(ed_s)+Fl;                                                    % Storing the local force vector into the global force vector
end
end