%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Thomas Wiltshire
% Date:   15/08/2018
% Description: Script to run example problem 3 (see documentation)
% Top tier of program. The user first specifies the problem
% type: either a problem with an anlytical solution, a problem with a
% non-analytical solution. This information is used to define the problem 
% geometry, boundary conditions to be imposed, type of mesh adaptivity to 
% use, number of adaptivity loops (if any), number of simulations to run 
% and material properties.
%
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
% 8. Post processing, plotting the error estimations etc.
% 
% Further information about the specific role of each function called can
% be found in the documentation.

%  Copyright (C) 2018 Thomas Wiltshire
%  $Revision: 1.0 $Date: 2018/08/15 17:09:20 $



% Test while loop----------------------------------------------------------
exit1 = 0;
while exit1==0
    % Problem definition --------------------------------------------------
    [node,edge,BC,delta2,delta1,loop_end,av_bc,sim_end,E,v]=...        % Problem definition for analytical problems
        nonanalytical_problem_3_generator;

    % Simulation while loop -----------------------------------------------
    for simulation = 1:sim_end
        DG=zeros(size(sim_end));L2=DG; Er=DG;ndof=DG;                      % Reset size of test data arrays
        
        % 1: Generating Mesh ----------------------------------------------
        [coord,etpl,etpl_face] = crackmesh(node,edge,BC);             % Generate mesh
        
        boundary_condition_check(etpl,etpl_face,coord);                    % Plot applied boundary conditions
        
        loop_count=0;  exit=0; nels=zeros(1,loop_end(simulation));         % Loop counter, adaptivity while loop exit condition, number of elements
        ed=[];                                                             % Empty variables for element dof topology and stiffness matrix
        fprintf('\n\nStart of simulation %d, delta 2=%.4f, delta 1=%.4f\n',...
            simulation,delta2(simulation),delta1(simulation));
        
        % 2: Start of adaptivity while loop -------------------------------
        while exit == 0
            loop_count=loop_count+1;                                                          % Number of adaptivity loops
            nels(loop_count)=sum(etpl.tree(:,1));                                             % Number of active elements
            fprintf('\n\n---------------------------------------------------------\n');
            fprintf('loop number:  %d | number of active elements: %d | max order:  %d\n',...
            loop_count,nels(loop_count),max(etpl.poly(:,2)));                                 % Loop counter
            etpl.ed_recalc(etpl.tree(:,1)==1)=1;
            
            % 3: Formulating the global stiffness matrix ------------------
            stiff_time=tic;fprintf('Stiffness calculation')
            [K,ndof(loop_count),ed] = DG_algorithm(etpl,etpl_face,coord,E,v);
            fprintf('\t%d s\n',toc(stiff_time));
            
            % 4: Force integration ----------------------------------------
            stiff_time=tic;fprintf('Force calculation')
            F_Dirichlet        = force_integration_Dirichlet(etpl,etpl_face,ed,coord,E,v);
            F_Dirichlet_Neuman = force_integration_Dirichlet_Neumann(etpl,etpl_face,ed,coord,E,v);
            F_Neumann          = force_integration_Neumann_nonanalytical(etpl,etpl_face,ed,coord);
            F_Body             = force_integration_vol(etpl,ed,coord,E,v);
            F=F_Body+F_Neumann+F_Dirichlet+F_Dirichlet_Neuman;

            fprintf(' \t    %d s\n',toc(stiff_time));
            
            % 5: Finding the solution -------------------------------------
            stiff_time=tic;fprintf('Solve calculation');
            u = force_solve(F,K,av_bc,etpl,coord,ed);
            fprintf(' \t    %d s\n',toc(stiff_time));
            
            %6: Error estimation-------------------------------------------
            
            % Error estimator ---------------------------------------------
            [El_Err,Er(loop_count)]=error_calc_ele_nonanalytical(u,etpl,ed,etpl_face,coord,E,v); 
            fprintf('Er error \t\t        %d\n',Er(loop_count))
            
            % Error in L2 norm --------------------------------------------
            [L2(loop_count)] = L2_norm_nonanalytical(coord,etpl,ed,u); 
            
            % Error in DG norm --------------------------------------------
            [DG(loop_count)] = DG_norm_calc_ele_nonanalytical(u,etpl,ed,etpl_face,coord,E,v); 
      
            
            % 7: hp-adaptivity --------------------------------------------
            if loop_count ~= loop_end(simulation)
                apt_time=tic;fprintf('Adapting mesh')
                [coord,etpl,etpl_face]=hp_adapt(El_Err,etpl,etpl_face,coord,delta2(simulation),delta1(simulation));
                fprintf(' \t\t    %d s\n',toc(apt_time));
            else
                exit = 1;
            end
        end
        
        % 8: Post-Processing-----------------------------------------------
        % Plot error estimators on the same axes---------------------------
        if loop_end(simulation)>1
            figure
            loglog(sqrt(ndof),Er);
            ndof_2=sqrt(ndof);
            s=sprintf('%d',-log(Er(end)/Er(end-1))/(log(ndof_2(end)/ndof_2(end-1))));
            text(mean(ndof_2),mean(Er),s);
            
            hold on;
            loglog(sqrt(ndof),DG);
            ndof_2=sqrt(ndof);
            s=sprintf('%d',-log(DG(end)/DG(end-1))/(log(ndof_2(end)/ndof_2(end-1))));  % Divide or minus
            text(mean(ndof_2),mean(DG),s);
            
            hold on;
            loglog(sqrt(ndof),L2);
            ndof_2=sqrt(ndof);
            s=sprintf('%d',-log(L2(end)/L2(end-1))/(log(ndof_2(end)/ndof_2(end-1))));
            text(mean(ndof_2),mean(L2),s);
            set(gca,'FontSize',16);
            
            legend('Error Estimate','DG norm error','L2 norm error');
            xlabel('ndof^{1/2}');
            ylabel('Error')
            
            disp_plotter(etpl,coord,ed,u,0.1);
            stressplotter_tri(E,v,coord,etpl,u,ed);
        
        end
        
    end

    
    exit1=1;
    
end