function testPrint(testFlag,simulation)  
% Function to print which test is being run
%
% Input(s):
% testFlag-          - Flag to indicate which test has been run
% simulation-        - Number of simulations flag
%
% Output(s):
%

%  Copyright (C) 2018 Thomas Wiltshire
%  $Revision: 1.0 $Date: 2018/08/10 17:09:20 $

    if testFlag==1 && simulation==1
        fprintf('Running the homogenous Dirichlet test problem with uniform h refinement\n');
    elseif testFlag==1 && simulation==2
        fprintf('Running the homogenous Dirichlet test problem with uniform p refinement\n');
    elseif testFlag==1 && simulation==3
        fprintf('Running the homogenous Dirichlet test problem with adaptive h refinement\n');    
    elseif testFlag==1 && simulation==4
        fprintf('Running the homogenous Dirichlet test problem with adaptive p refinement\n'); 
    elseif testFlag==1 && simulation==5
        fprintf('Running the homogenous Dirichlet test problem with adaptive hp refinement\n'); 
    elseif testFlag==2
        fprintf('Running the inhomogenous Dirichlet test problem\n');
    elseif testFlag==3
        fprintf('Running the inhomogenous Neumann test problem\n');
    elseif testFlag==4
        fprintf('Running the inhomogenous mixed test problem\n');
    end
end