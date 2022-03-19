function testProblemCheck(simulation,sim_end,testFlag,Er,DG,L2,ndof,nels)
% Function to display results of test problems
%
% Input(s):
% simulation-        - Number of simulations flag
% sim_end            - Number of simulations to run
% testFlag           - Flag to indicate which test has been run
%
% Ouput(s);
% 

%  Copyright (C) 2018 Thomas Wiltshire
%  $Revision: 1.0 $Date: 2018/08/10 17:09:20 $

fprintf('\n\n---------------------------------------------------------\n');
if testFlag==1 && simulation ==1
    fprintf('\nResults of Test 1 with pure h refinement, u=v=sin(pi*x)sin(pi*y)\nDG Norm:        %.16f  %.16f  %.16f\nL2 Norm:        %.16f  %.16f  %.16f\nError Estimate: %.16f  %.16f  %.16f\nNdof:           %d %d %d\nNo. Elements:   %d  %d   %d\n',DG,L2,Er,ndof,nels);  
elseif testFlag==1 && simulation ==2
    fprintf('\nResults of Test 1 with pure p refinement, u=v=sin(pi*x)sin(pi*y)\nDG Norm:        %.16f  %.16f  %.16f\nL2 Norm:        %.16f  %.16f  %.16f\nError Estimate: %.16f  %.16f  %.16f\nNdof:           %d %d %d\nNo. Elements:   %d  %d  %d\n',DG,L2,Er,ndof,nels);
elseif testFlag==1 && simulation ==3
    fprintf('\nResults of Test 1 with adaptive h refinement, u=v=sin(pi*x)sin(pi*y)\nDG Norm:        %.16f  %.16f  %.16f\nL2 Norm:        %.16f  %.16f  %.16f\nError Estimate: %.16f  %.16f  %.16f\nNdof:           %d %d %d\nNo. Elements:   %d  %d  %d\n',DG,L2,Er,ndof,nels);
elseif testFlag==1 && simulation ==4
    fprintf('\nResults of Test 1 with adaptive p refinement, u=v=sin(pi*x)sin(pi*y)\nDG Norm:        %.16f  %.16f  %.16f\nL2 Norm:        %.16f  %.16f  %.16f\nError Estimate: %.16f  %.16f  %.16f\nNdof:           %d %d %d\nNo. Elements:   %d  %d  %d\n',DG,L2,Er,ndof,nels);
elseif testFlag==1 && simulation ==5
    fprintf('\nResults of Test 1 with h p refinement, u=v=sin(pi*x)sin(pi*y)\nDG Norm:        %.16f  %.16f  %.16f\nL2 Norm:        %.16f  %.16f  %.16f\nError Estimate: %.16f  %.16f  %.16f\nNdof:           %d %d %d\nNo. Elements:   %d  %d  %d\n',DG,L2,Er,ndof,nels);
    
elseif testFlag==2
    fprintf('\nResults of Test 2, inhomogenous Dirichlet boundary conditions, u=sin(pi*x)sin(pi*y), v=cos(pi*x)sin(pi*y)\nDG Norm:        %.16f\nL2 Norm:        %.16f\nError Estimate: %.16f\nNdof:           %d\nNo.Elements:    %d\n',DG,L2,Er,ndof,nels);
    
elseif testFlag==3
    fprintf('\nResults of Test 3, inhomogenous Neumann boundary conditions, u=sin(2pi*x)sin(2pi*y), v=cos(2pi*x)sin(2pi*y)\nDG Norm:        %.16f\nL2 Norm:        %.16f\nError Estimate: %.16f\nNdof:           %d\nNo.Elements:    %d\n',DG,L2,Er,ndof,nels);

elseif testFlag==4
    fprintf('\nResults of Test 4, mixed boundary conditions, u=x, v=y\nDG Norm:        %.16f\nL2 Norm:        %.16f\nError Estimate: %.16f\nNdof:           %d\nNo.Elements:    %d\n',DG,L2,Er,ndof,nels);
end

if simulation==sim_end
    fprintf('End of simulation\n')
elseif sim_end>1
    fprintf('End of simulation %d.\n', simulation)
    
end

