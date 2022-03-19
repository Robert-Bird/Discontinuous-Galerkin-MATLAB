function [rot]=rot_ave_BC(etpl,coord,ed)
% Compute rotational components of the force.
% Returns Lagrangian 'rotational force' component to force a system to act
% in a deterministic manner
%
% Input(s):
% etpl       - Element tolopogy struture 
% coord      - Element coordinates
% ed         - Degree of freedom steering matrix
%
% Ouput(s):
% rot - Lagrangian rotational 'force' rotationa about the z direction
%
%  See also UV_AVE_BC.

% Accounts for varying element polynomial order
p_max=max(etpl.poly(etpl.tree(:,1)==1,2)); nD = 2;
[ngp,vwp,dNr,~]=dNr_vol_gauss(p_max);
rot=zeros(max(ed(:)),1);
nels = size(etpl.mat,1);
nD=2;

for nel = 1:nels
    if etpl.tree(nel,1)==1
        loc_p=etpl.poly(nel,2);
        nov_p=nov_calc(loc_p);
        rot2=zeros(nov_p*2,1);
        dNxt=zeros(nov_p*2,1);
        JT=dNr(:,1:3)*coord(etpl.mat(nel,:),:);
        dNrp=choose_dNr(dNr,loc_p,p_max);
        for gp = 1:ngp
            indx=2*(gp-1)+(1:nD);
            dNx=JT(indx,:)\dNrp(indx,:);
            dNxt(1:2:end)=dNx(2,:);
            dNxt(2:2:end)=-dNx(1,:);
            detJ=det(JT(indx,:));
            rot2=rot2+dNxt*detJ*vwp(gp);
        end
        rot(ed(nel,1:(nov_p*2)),:)=rot2;
    end
end
