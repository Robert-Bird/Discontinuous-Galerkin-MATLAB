function F=force_integration_Dirichlet(etpl,etpl_face,ed,coord,E,v)
% Dirichlet boundary condition.
%
% Input(s):
% etpl       - Element tolopogy struture (see seed_mesh.m)
% etpl_face  - Element face tolopogy struture (see seed_mesh.m)
% ed         - Degrees of freedom steering matrix
% pr         - Pressure
% coord      - Element coordinates (see seed_mesh.m)
% Dirichlet  - Dirichlet BC function handle
% F_h        - Dirichlet BC function handle
% E          - Young's Modulus
% v          - Poisson's ratio
%
% Ouput(s):
% F - Force vector from Dirichlet boundary condition

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

F=zeros(max(max(ed)),1);                                                   % Assigning variable space for the force vector
etpl_face = etpl_face(etpl_face(:,end)==-2,:);                             % Face information for external elements
p_max=max(etpl.poly(:,2)); nD = 2;                                         % Maximum polynomial order for the external faces
[swp,dNr,Nr,ngp]=dNr_surf_gauss(p_max);                                    % Shape functions, derivatives and guass points
n_faces = size(etpl_face,1);                                               % Number of external faces
nx   = etpl_face(:,5);     ny = etpl_face(:,6);                            % Outward normal componenets and lenght of element face.
en   = etpl_face(:,1);   f_p  = etpl_face(:,3);                            % Positive element number and faces.
De=d_mat(E,v);                                                             % Hookian stiffness matrix
for i = 1:n_faces                                                          % Integral loop
    pos_loc_p=etpl.poly(en(i),2);                                          % Positive element polynomial order
    nov_p=nov_calc(pos_loc_p)*2;                                           % Positive element number of variables
    Fl = zeros(nov_p,1); 
    p_f=f_p(i);                                                            % Positive element face number
    Bp=zeros(3,nov_p);Npm=zeros(2,nov_p);
    nm = [nx(i) 0 ny(i);0 ny(i) nx(i)];                                    % Outward normal matrix
    h=face_size_cal(p_f,coord(etpl.mat(en(i),:),:));                       % Length of the face
    for gp = 1:ngp
        indx_dNr = ((gp-1)*nD)+((p_f-1)*nD*ngp)+(1:nD);                    % Positive shape function derivative index
        indx_sNr = gp+((p_f-1)*ngp);                                       % Positive shape function index
        co_pve   = coord(etpl.mat(en(i),:),:);                             % Coordinates of positive element
        Jp       = dNr(indx_dNr,1:3)*co_pve;                               % Jacobian of positive element
        [pen,~]      = pen_func(co_pve,co_pve,h,pos_loc_p,pos_loc_p,E,v,nx(i),ny(i));	% Penalty term
        dNrp     = choose_dNr(dNr,pos_loc_p,p_max); dNx_p=Jp\dNrp(indx_dNr,:);	% Global shape function derivatives
        Nrp      = choose_Nr(Nr,pos_loc_p,p_max); Np = Nrp(indx_sNr,:);    % Local shapen functions
        Bp([1 3],1:2:end) = dNx_p; Bp([3 2],2:2:end)=dNx_p;                % Shape function derivative matrix
        Npm(1,1:2:end)    = Np; Npm(2,2:2:end) = Np;                       % Shape function matrix positive face
        W                 = swp(gp)*h/2;                                   % Integral weight
        X=(Np(1,1:3)*coord(etpl.mat(en(i),:),1));                       % Global x-position
        Y=(Np(1,1:3)*coord(etpl.mat(en(i),:),2));                       % Global y-position
        u=ImposedDisplacements(X,Y);                                                  %F_h(nx(i),ny(i),X,Y,pr,E,v,U_func);    % Displacement vector at global position (x,y)
        Fl = Fl - ((Bp'*De*nm'))*W*u +  pen*(Npm')*W*u;            % Local stiffness calculation	% Local stiffness calculation
    end
    ed_s=ed(en(i),1:sum(ed(en(i),:)~=0));                                  % Local degrees of freedom
    F(ed_s)=F(ed_s)+Fl;                                                    % Storing the local force vector into the global force vector
end


