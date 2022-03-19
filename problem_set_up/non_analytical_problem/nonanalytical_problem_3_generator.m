function [node,edge,BC,d_2,d_1,loop_end,av_bc,sim_end,E,nu]=nonanalytical_problem_3_generator()
% Problem set up function
%
% Ouput(s);
% node-     - Vertices of the desired domain (for mesh generation)
% edge-     - Topology matrix of the desired domain (for mesh generation)
% BC-       - Boundary conditions to be  imposed:
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
% loop_end- - Number of adaptivity loops
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
node=[0 -1.5;                                                              % Node 1
      1.5 -1.5;                                                            % Node 2
      1.5 1.5;                                                             % Node 3
      0 1.5;                                                               % Node 4
      0 eps();                                                             % Node 5: Set crack lip nodes as almost coincident
      0.5 0;                                                               % Node 6
      0 -eps()                                                             % Node 7: Set crack lip nodes as almost coincident
      ];                                                                   
  
edge=[
    1 2;                                                                   % Edge 1
    2 3;                                                                   % Edge 2
    3 4;                                                                   % Edge 3
    4 5;                                                                   % Edge 4
    5 6;                                                                   % Edge 5
    6 7;                                                                   % Edge 6
    7 1];                                                                  % Edge 7

BC=[edge,3*ones(size(edge,1),1)];                                          % Define default boundary conditions
BC(1,3)=-3;
BC(3,3)=2;
BC(2,3)=4;
E = 2.5; nu = 0.25;                                                        % Young's Modulus (Pa); Poisson's Ratio

sim_end=3;                                                                 % Number of different simulations to run (for example with different delta values)
loop_end=[15 15 5];                                                        % Number of adaptivity loops

% Strategy 1 coefficients
d_2=[0.3 0.2 0];                                                           % hp adaptive, h adaptive, pure h
d_1=[0.07 0.2 0];                                                           

syms X Y
u=0*X;                                                                     % Displacements imposed
v=0*Y;

traction= [0*X; -1];

fx=0*X;fy=0*X;

if isempty(d_2)|| isempty(d_1)
    error('No input for strategy 1 coefficients')
elseif max(size(loop_end))~=sim_end
    warning('Number of adaptivity loops specified does not match number of simulations to be performed')
end

if isempty(BC)==0
    if all(BC(:,3)==-3)
        av_bc=1;
    else
        av_bc=0;
    end
else
    av_bc=0;
end
De=d_mat(E,nu);

ImposedDisplacements(X,Y)=[u;v];
matlabFunction(ImposedDisplacements,'File','problem_set_up/ImposedDisplacements');        % Create displacement function using Matlab symbolic toolbox

% Stress

traction(X,Y)=traction;
matlabFunction(traction,'File','problem_set_up/tractionFunc');                            % Create traction function using Matlab symbolic toolbox

% Body force
BodyForce(X,Y)=[fx;fy];
matlabFunction(BodyForce,'File','problem_set_up/BodyForce');                              % Create body force function using Matlab symbolic toolbox
end