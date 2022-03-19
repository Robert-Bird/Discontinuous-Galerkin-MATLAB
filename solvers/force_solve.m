function u=force_solve(f,K,av_bc,etpl,coord,ed)
% Solve the linear system.
% Check if a node exists.
% Returns the displacement solution.
%
% Input(s):
% f     - Force vector of the mesh
% K     - SIPG stiffness matrix
% av_bc - Average BC flag
% u_bc  - u displacement lagrangian boundary conditions
% v_bc  - v displacement lagrangian boundary conditions
% etpl  - Element topology structure
% coord - Element nodal coordinates
% ed    - Element degree of freedom matrix
%
% Ouput(s):
% u - displacement from solving linear system of equations

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 

if av_bc == 1
    [u_bc,v_bc] = UV_ave_BC(etpl,coord,ed);                                % Displacement average boundary conditions
    [rot] = rot_ave_BC(etpl,coord,ed);                                    % Rotation average boundary condition
    f=[f',0 0 0]'; K=[K,u_bc,v_bc,rot]; K=[K;[u_bc',0,0,0];[v_bc',0,0,0];[rot',0,0,0]];
    u=K\f;
else
    K=(K'+K)./2;
    u=K\f;
end
