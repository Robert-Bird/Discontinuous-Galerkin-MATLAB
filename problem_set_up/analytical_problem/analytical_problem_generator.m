function [node,edge,BC,d_2,d_1,loop_end,av_bc,sim_end,E,nu]=analytical_problem_generator()
% Problem set up function
% Analytical problem generation. Using Matlab's symbolic math toolbox,
% generates functions for displacement, stress and body force across
% the domain. The stress expression is generated by differentiating the
% displacement with respect to cartesian coordinates and multiplying by
% the elastic stiffness matrix, see d_mat function. To find the body force
% associated with this displacement, the stress is then differentiated
% w.r.t the cartesian corrdinates.
%
% Problem set up script for problems with an analytical 
% solution. The user should first specify the boundary conditions, 
% problem geometry and material properties. To specify the problem 
% geometry, an input of the [x,y] coordinates of the vertices of the
% domain is required,for example a unit square would have 
% nodes=[0,0;1,0;1,1,0,1].
% 
% The default boundary condition is homogenous Neumann, indicated by a 3 in
% the third column of BC. To impose a boundary condition simply choose the
% appropriate row in the BC array; note that the face numbering is
% anticlockwise. To check the correct problem geometry has been defined 
% and the correct BCs have been imposed, a plot of the geometry
% has been included in the code with external faces colour coded according
% to the boundary condition on that face.
% 
% An input is then required for the desired method of mesh adaptivity, 
% using the values of delta 2 and delta 1. Should it be desired that
% multiple simulations with different methods of mesh adaptivity are 
% to be performed, arrays of delta values can be specified for example 
% d_2=[1 0.5 0.3] d_1=[1 0.45 0.2]. Note that delta 2 must always be
% greater than delta 1.
%
% Finally, as this script is for problems where the displacement solution
% is known across the domain, an expression for displacement must be
% created using the Matlab symbolic toolbox. This is then passed on to the
% function analytical_problem_generator where an expression for the
% analytical body force associated with the analytical displacement
% solution is generated, i.e the RHS {f} of the weak form of equilibrium
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

node = [                                                                   % list of xy "node" coordinates (see Mesh2D documentation)
    0,0
    2,0
    0,2
    ] ;

edge = [                                                                   % list of "edges" between nodes (see Mesh2D documentation)
    1, 2
    2, 3
    3, 1
    ] ;

BC=[edge,3*ones(size(edge,1),1)];                                          % Default boundary condition

BC(:,3)=2;
BC(2, 3)=-2;

E = 5/2; nu = 0.25;                                                        % Young's Modulus (Pa); Poisson's Ratio

sim_end=2;                                                                 % Number of different simulations to run (for example with different delta values)
loop_end=[2 2];                                                                % Number of adaptivity loops

% Strategy 1 coefficients
d_2= [0.7 0];                                                              % delta 2 
d_1=[0.07 0];                                                              % delta 1


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