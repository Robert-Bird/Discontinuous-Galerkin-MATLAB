function F=force_integration_Neumann_nonanalytical(etpl,etpl_face,ed,coord)
% Force integration for Neumann faces.
% Currently only a force integration for the Neumann boundary conditions
% this function needs to be expanded to included all force types
%
% Input(s):
% etpl       - Element tolopogy struture (see seed_mesh.m)
% etpl_face  - Element face tolopogy struture (see seed_mesh.m)
% ed         - Degrees of freedom steering matrix
% pr         - Pressure
% coord      - Element coordinates (see seed_mesh.m)
%
% Ouput(s):
% F - Force vector from Neumann boundary condition

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

F=zeros(max(max(ed)),1);                                                   % Assigning variable space for the force vector
p_max=max(etpl.poly(:,2));                                                 % Maximum polynomial order of elements which have a Newmann boundary condition
etpl_face=etpl_face(etpl_face(:,end)==-3,:);                               % Newmann boundary conditions that the Newmann boundary condition exists
[swp,~,Nr,ngp]=dNr_surf_gauss(p_max);                                      % Generating shape functions and Gauss point weights
if size(etpl_face,1)>0
    en = etpl_face(:,1);                                                   % Element number
    fn = etpl_face(:,3);                                                   % Face number
    nx = etpl_face(:,5);                                                   % normal x
    ny = etpl_face(:,6);                                                   % normal y
    for i = 1:size(etpl_face,1)
        h=face_size_cal(fn(i),coord(etpl.mat(en(i),:),:));
        Fl  = zeros(sum(ed(en(i),:)~=0),1);                                % Local force component
        Npm  = zeros(2,sum(ed(en(i),:)~=0));                               % Assigning memory for shape function matrix
        loc_p=etpl.poly(en(i),2);                                          % positive element polynomial
        for gp = 1:ngp
            indx_sNr = gp+((fn(i)-1)*ngp);                                 % Shape function index
            Np  = choose_Nr(Nr,loc_p,p_max);                               % Choosing shape functions
            Npm(1,1:2:end) = Np(indx_sNr,:);
            Npm(2,2:2:end) = Np(indx_sNr,:);
            W = swp(gp)*h/2;                                               % Combining Gauss point with Jacobian
            X=(Npm(1,1:2:6)*coord(etpl.mat(en(i),:),1));                   % Global x-position
            Y=(Npm(1,1:2:6)*coord(etpl.mat(en(i),:),2));                   % Global y-position
            Fl=Fl+Npm'*tractionFunc(X,Y).*W;   % Summing the force vector result
            
        end
        ed_s=ed(en(i),1:sum(ed(en(i),:)~=0));                              % Local degrees of freedom
        F(ed_s)=F(ed_s)+Fl;                                                % Storing the local force vector into the global force vector
    end
    
end
end
