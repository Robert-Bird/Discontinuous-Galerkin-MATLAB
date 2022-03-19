function [u,v]=UV_ave_BC(etpl,coord,ed)
% Compute components of the force.
% Returns Lagrangian force components to force a system to act 
% in a deterministic manner.
%
% Input(s):
% etpl       - Element tolopogy struture (see seed_mesh.m)
% coord      - Element coordinates (see seed_mesh.m)
% ed         - Degree of freedom steering matrix
%
% Ouput(s):
% u - Lagrangian 'force' for u displacement
% v - Lagrangian 'force' for v displacement
%
%  See also ROT_AVE_BC.

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $
p_max=max(etpl.poly(etpl.tree(:,1)==1,2)); nD = 2;
[ngp,vwp,dNr,Nr]=dNr_vol_gauss(p_max);
u=zeros(max(ed(:)),1);v=u;
nels = size(etpl.mat,1);
nD=2;
for nel = 1:nels
    if etpl.tree(nel,1)==1
        loc_p=etpl.poly(nel,2);                    
        nov_p=nov_calc(loc_p);
        u2=zeros(nov_p,1);
        JT=dNr(:,1:3)*coord(etpl.mat(nel,:),:);
        Nrp=choose_Nr(Nr,loc_p,p_max);
        for gp = 1:ngp
            indx=2*(gp-1)+(1:nD);
            Nx=Nrp(gp,:)';
            detJ=det(JT(indx,:));
            u2=u2+Nx*detJ*vwp(gp);
        end
        u(ed(nel,1:2:nov_p*2),:)=u2;
        v(ed(nel,2:2:nov_p*2),:)=u2;
    end
end