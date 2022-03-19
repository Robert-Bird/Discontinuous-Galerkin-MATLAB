function stressplotter_tri(E,v,coord,etpl,u,ed)
% Creates a surface plot of the von mises stress for a 2D problem
%
% Input(s):
% E          - Youngs modulus
% v          - Poisson's ratio
% coord      - Element coordinates 
% etpl       - Element tolopogy struture 
% u          - Displacement solution created by the SIPG algorithm
% ed         - Degrees of freedom steering matrix

%  Copyright (C) 2018 Robert Bird
%  $Revision: 1.0 $Date: 2018/08/21 17:09:20 $

nD = 2;                                                                    % Number of dimensions in the problem
p_max=max(etpl.poly(etpl.tree(:,1)==1,2));                                 % Maximum polynomial order of the current mesh
[~,~,dNr,Nr]=dNr_vol_gauss_plotter(p_max);                                 % Volume shape functions, derivatives and guass points                               
[~,dNrs,Nrs,~]=dNr_surf_gauss(p_max);                                      % Surface shape functions, derivatives and guass points
dNr= [dNr;dNrs];                                                           % Putting all shape function derivatives into the same matrix
Nr = [Nr;Nrs];                                                             % Putting all shape functions into the same matrix
ngp=size(Nr,1);                                                            % Number of Gauss points to loop through
De=d_mat(E,v);                                                             % Hookian stiffness matrix
nels = size(etpl.mat,1);                                                   % Number of elements in the mesh
sigma=zeros(nels,3);                                                       % matrix of [von mises stress,x-coordinate,y-coordinate]
tri=zeros(10*nels*ngp,3);                                                  % Triangulation matrix for total mesh
c=0;                                                                       % Counter for tri_temp storage in tri (line 56)                                                                 
stress_count=0;                                                            % Counter for storage of von-mises stress and coordinates into sigma
nel_count=0;                                                               % Counter for the number of elements over which the stress is calculated 
for nel = 1:nels                                                           % Integral loop over all elements
    if etpl.tree(nel,1)==1                                                 % Checking if element was used in the current mesh
        nel_count=nel_count+1;                                             % Increasing the element number counter                     
        loc_p=etpl.poly(nel,2);                                            % Local polynomial order
        nov=nov_calc(loc_p)*2;                                             % Number of variables
        JT=dNr(:,1:3)*coord(etpl.mat(nel,:),:);                            % Jocabian
        dNrt=choose_dNr(dNr,loc_p,p_max);                                  % Local shape function derivatives                             
        B=zeros(3,nov);                                                    % Assigning space for derivative matrix
        for gp = 1:ngp                                                     % Gauss point loop
            indx=2*(gp-1)+(1:nD);                                          % Jacobian and shape function derivative index
            dNx=JT(indx,:)\dNrt(indx,:);                                   % Global shape function derivative
            B([1 3],1:2:end)=dNx;                                          % Derivative matrix
            B([3 2],2:2:end)=dNx;                                          %
            strain=B*u(ed(nel,ed(nel,:)~=0));                              % Strain Calculation
            s=De*strain;                                                   % Cauchy stress
            sm=[s(1) s(3);s(3) s(2)];                                      % Forming a Cauchy stress matrix
            s_max=max(eig(sm))-min(eig(sm));                               % Von-Mises stress calculation
            stress_count=stress_count+1;                                   % Increasing the stress storage counter
            X = Nr(gp,1:3)*(coord(etpl.mat(nel,1:3),1));                   % X - coordinate
            Y = Nr(gp,1:3)*(coord(etpl.mat(nel,1:3),2));                   % Y - coordinate
            sigma(stress_count,:) =  [s_max,X,Y];                          % Storage of Von-Mises stress and global coordinates
        end
        tri_temp = delaunay(sigma(ngp*(nel_count-1)+(1:ngp),2),sigma(ngp*(nel_count-1)+(1:ngp),3)); % Determining triangulation for element nel
        tri(c+(1:size(tri_temp,1)),:) = tri_temp+max(tri(:));              % Storing the triangulation
        c=c+size(tri_temp,1);                                              % Keeping track of the number of rows used in tri
    end
end
tri((c+1):end,:)=[];                                                       % Removing excess lines of tri
figure                                                                     % Creating a new figure
trisurf(tri,sigma(:,2),sigma(:,3),sigma(:,1),'EdgeColor','none','LineStyle','none','FaceLighting','phong'); % Plotting stress
xlabel('x coordinate (m)');                                                    % Labels
ylabel('y coordinate (m)');                                                    %
zlabel('Stress (Pa)');                                                          %


