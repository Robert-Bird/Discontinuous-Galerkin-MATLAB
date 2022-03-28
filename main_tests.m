% Main script for running test problems.
% The program then performs the number of simulations specified. A
% simulation consists of:
% 1. Generating the mesh using Mesh2D software
% 2. Entering an adaptivity while loop
% 3. Generating a discontinuous Galerkin stiffness matrix (see
%    documentation)
% 4. Calculating force vectors to balance the weak form of the linear
%    elastic problem (see documentation)
% 5. Calculating an FE approximation of the solution
% 6. Calculating error estimations (see documentation)
% 7. Refining the mesh according to an adaptivity algorithm and user
%    specified values of delta 1 and delta 2
%
% Further information about the specific role of each function called can
% be found in the documentation.

%  Copyright (C) 2018 Thomas Wiltshire
%  $Revision: 1.0 $Date: 2018/08/15 17:09:20 $


% Defines the first test that is performed

% Counter of passes tests 
%(if testFlag=1, complete success is a values of 5)
%(if testFlag=2, complete success is a values of 1)
%(if testFlag=3, complete success is a values of 1)
%(if testFlag=4, complete success is a values of 1)



passcheck = 0;

% Defines the first test that is performed
global testFlag


% Problem definition ------------------------------------------------------
% Test problem definitions
[sim_end,loop_end,testFlag,E,v,delta2,delta1,av_bc]=...
    testProblem(testFlag);

% Simulation while loop ---------------------------------------------------
for simulation = 1:sim_end
    DG=zeros(size(sim_end));L2=DG; Er=DG;ndof=DG;                      % Reset size of test data arrays
    
    % 1: Generating Mesh --------------------------------------------------
    [coord,etpl,etpl_face] = testmesh(testFlag);                       % Generate mesh
    loop_count=0;  exit=0; nels=zeros(1,loop_end(simulation));         % Loop counter, adaptivity while loop exit condition, number of elements
    ed=[];                                                             % Empty variables for element dof topology and stiffness matrix
    testPrint(testFlag,simulation);                                    % State which test problem is being performed (if test_switch==1)
    fprintf('\n\nStart of simulation %d, delta 2=%.4f, delta 1=%.4f\n',...
        simulation,delta2(simulation),delta1(simulation));
    
    % 2: Start of adaptivity while loop -----------------------------------
    while exit == 0
        loop_count=loop_count+1;                                                          % Number of adaptivity loops
        nels(loop_count)=sum(etpl.tree(:,1));                                             % Number of active elements
        fprintf('\n\n---------------------------------------------------------\n');
        fprintf('loop number:  %d | number of active elements: %d | max order:  %d\n',...
            loop_count,nels(loop_count),max(etpl.poly(:,2)));                                 % Loop counter
        etpl.ed_recalc(etpl.tree(:,1)==1)=1;
        
        % 3: Formulating the global stiffness matrix ----------------------
        stiff_time=tic;fprintf('Stiffness calculation')
        [K,ndof(loop_count),ed] = DG_algorithm(etpl,etpl_face,coord,E,v);
        fprintf('\t%d s\n',toc(stiff_time));
        
        % 4: Force integration --------------------------------------------
        stiff_time=tic;fprintf('Force calculation')
        F_Dirichlet        = force_integration_Dirichlet(etpl,etpl_face,ed,coord,E,v);
        F_Dirichlet_Neuman = force_integration_Dirichlet_Neumann(etpl,etpl_face,ed,coord,E,v);
        F_Neumann          = force_integration_Neumann_analytical(etpl,etpl_face,ed,coord);
        F_Body             = force_integration_vol(etpl,ed,coord,E,v);
        F=F_Body+F_Neumann+F_Dirichlet+F_Dirichlet_Neuman;
        
        fprintf(' \t    %d s\n',toc(stiff_time));
        
        % 5: Finding the solution -----------------------------------------
        stiff_time=tic;fprintf('Solve calculation');
        u = force_solve(F,K,av_bc,etpl,coord,ed);
        fprintf(' \t    %d s\n',toc(stiff_time));
        
        %6: Error estimation-----------------------------------------------
        
        % Error estimator -------------------------------------------------
        [El_Err,Er(loop_count)]=error_calc_ele_analytical(u,etpl,ed,etpl_face,coord,E,v);
        fprintf('Er error \t\t        %d\n',Er(loop_count))
        
        % Error in L2 norm ------------------------------------------------
        [L2(loop_count)] = L2_norm_analytical(coord,etpl,ed,u);
        
        % Error in DG norm ------------------------------------------------
        [DG(loop_count)] = DG_norm_calc_ele_analytical(u,etpl,ed,etpl_face,coord,E,v);
        
        
        % 7: hp-adaptivity ------------------------------------------------
        if loop_count ~= loop_end(simulation)
            apt_time=tic;fprintf('Adapting mesh')
            [coord,etpl,etpl_face]=hp_adapt(El_Err,etpl,etpl_face,coord,delta2(simulation),delta1(simulation));
            fprintf(' \t\t    %d s\n',toc(apt_time));
        else
            exit = 1;
        end
    end
    
    
    [passcheckNew]=testCheck(simulation,sim_end,testFlag,Er,DG,L2,ndof,nels,passcheck);
    passcheck=passcheckNew;
end



