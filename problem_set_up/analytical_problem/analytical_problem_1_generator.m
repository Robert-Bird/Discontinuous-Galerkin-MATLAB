function [node,edge,BC,d_2,d_1,loop_end,av_bc,sim_end,E,nu]=analytical_problem_1_generator()
% Problem set up function
% Analytical problem generation. Using Matlab's symbolic math toolbox,
% generates functions for displacement, stress and body force across
% the domain. The stress expression is generated by differentiating the
% displacement with respect to cartesian coordinates and multiplying by
% the elastic stiffness matrix, see d_mat function. To find the body force
% associated with this displacement, the stress is then differentiated
% w.r.t the cartesian corrdinates.
%
% Ouput(s);
% node     - Vertices of the desired domain (for mesh generation)
% edge     - Topology matrix of the desired domain (for mesh generation)
% BC       - Boundary conditions to be  imposed:
%
% Boundary Condition               Flag in BC
% Homogeneous Dirichlet             2
% Prescribed  Dirichlet            -2
% Homogeneous Neumann               3
% Prescribed Neumann               -3
% Homogenous Dirichlet/Neumann      4
% Prescribed Dirichlet/Neumann     -4
%
% d_max-    - Delta 2 coefficient: Used to determine which kind of mesh
%             refinement regime to use
% d_min-    - Delta 1 coefficient: Used to determine which kind of mesh
%             refinement regime to use
%
% hp refinement strategy   d_min        d_max           Comment
% Pure h                   0            0               -------
% Pure p                   0            1               -------
% Adaptive h               0=<d_min=<1  0=<d_max=<1     d_min=d_max
% Adaptive p               0=<d_min=<1  1               -------
% hp Adaptive              0=<d_min=<1  0=<d_max=<1     d_max>d_min
% No adaptivity            1            1               -------
%
% loop_end- - Number of adaptivity loops array, size of which must equal
%             sim_end
% av_bc-    - Average boundary condition flag; turn on when all
%             boundary conditions are inhomogenous Neumann
% sim_end-  - Number of different simulations to run (default=1)
% testFlag- - Flag to indicate which test has been run (default=0)
%
% Test Problem                      testFlag Value
% hp Refinement                     1
% Inhomogenous Dirichlet            2
% Inhomogenous Neumann              3
% Mixed Dirichlet/ Neumann          4
% testFlag-
% E-        - Young's Modulus (GPa)
% nu-       - Poisson's Ratio

%  Copyright (C) 2018 Thomas Wiltshire
%  $Revision: 1.0 $Date: 2018/08/10 17:09:20 $

% list of xy "node" coordinates 
node = [                                                                   
    0,0                                                                    % Node 1
    1,0                                                                    % Node 2
    1,1                                                                    % Node 3
    0,1                                                                    % Node 4
    ] ;

% list of "edges" between nodes
edge = [                                                                   
    1, 2                                                                   % Edge 1
    2, 3                                                                   % Edge 2
    3, 4                                                                   % Edge 3
    4, 1                                                                   % Edge 4
    ] ;

BC=[edge,3*ones(size(edge,1),1)];                                          % Default boundary condition inhomogenous Neumann
BC(:,3)=2;                                                                 % Apply homogenous Dirichlet across all faces

E = 5/2; nu = 0.25;                                                        % Young's Modulus (Pa); Poisson's Ratio

sim_end=3;                                                                 % Number of different simulations to run (for example with different delta values)
loop_end=[11 4 3];                                                         % Number of adaptivity loops

% Strategy 1 coefficients
d_2= [0.7 0.07 0];                                                         % hp adaptive, h adaptive, uniform h 
d_1=[0.07 0.07 0] ;                                                         

syms X Y
u=sin(2*pi*X)*sin(2*pi*Y);                                                 % Displacement solution
v=sin(2*pi*X)*sin(2*pi*Y);


if isempty(d_2)|| isempty(d_1)
    error('No input for strategy 1 coefficients')
elseif max(size(loop_end))~=sim_end
    warning('Number of adaptivity loops specified does not match number of simulations to be performed')
end

if all(BC(:,3)==-3)
    av_bc=1;
else
    av_bc=0;
end
De=d_mat(E,nu);

ImposedDisplacements(X,Y)=[u;v];
matlabFunction(ImposedDisplacements,'File','problem_set_up/ImposedDisplacements');        % Create Matlab displacement function from symbolic variables

% Diff of displacement
du_dx=diff(u,X);
du_dy=diff(u,Y);

dv_dx=diff(v,X);
dv_dy=diff(v,Y);

% Strain
e = [du_dx; dv_dy; (du_dy+dv_dx)];

% Stress
s = simplify(De*e);
BodyStress(X,Y)=s;
matlabFunction(BodyStress,'File','problem_set_up/BodyStress');                            % Create Matlab stress function from symbolic variables

% Body force
fx = - simplify(diff(s(1),X)+diff(s(3),Y));
fy = - simplify(diff(s(2),Y)+diff(s(3),X));
BodyForce(X,Y)=[fx;fy];
matlabFunction(BodyForce,'File','problem_set_up/BodyForce');                              % Create Matlab body force function from symbolic variables
end