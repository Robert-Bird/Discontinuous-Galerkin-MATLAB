function [sim_end,loop_end,testFlag,E,nu,d_max,d_min_poly,av_bc]=testProblem(testFlag)
% Test set up function.
%
% Input(s):
% test_flag-     - Indicates which test problem to run
%
% Ouput(s);
% BC-            - Boundary conditions
% sim_end-       - Number of different simulations to run
% loop_end-      - Number of adaptivity loops
% node-          - Vertices of the desired domain (for mesh generation)
% edge-          - Topology matrix of the desired domain (for mesh generation)
% E-             - Young's Modulus (GPa)
% nu-            - Poisson's Ratio
% Test Problem                      testFlag Value
% hp Refinement                     1
% Inhomogenous Dirichlet            2
% Inhomogenous Neumann              3
% Mixed Dirichlet/ Neumann          4

%  Copyright (C) 2018 Thomas Wiltshire
%  $Revision: 1.0 $Date: 2018/08/10 17:09:20 $

E = 1; nu = 0.3;                                                           % Young's Modulus; Poisson's Ratio
De=d_mat(E,nu);                                                            % Hookian stiffness matrix

syms X Y;

if testFlag==1
    u=sin(pi*X)*sin(pi*Y);
    v=u;
    av_bc = 0;
    sim_end=5;
    loop_end=[3 3 3 3 3];
    d_max=[0 1 0.3 1 0.7]; d_min_poly=[0 0 0.3 0.3 0.3];
elseif testFlag==2
    u=sin(pi*X)*sin(pi*Y);
    v=cos(pi*X)*sin(pi*Y);
    av_bc = 0;
    sim_end=1;
    loop_end=1;
    d_max=1; d_min_poly=1;
elseif testFlag==3
    u=sin(2*pi*X)*sin(2*pi*Y);
    v=cos(2*pi*X)*sin(2*pi*Y);
    av_bc = 1;
    sim_end=1;
    loop_end=1;
    d_max=1; d_min_poly=1;
elseif testFlag==4
    u=X;
    v=Y;
    av_bc = 0;
    sim_end=1;
    loop_end=1;
    d_max=1; d_min_poly=1;
end


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
