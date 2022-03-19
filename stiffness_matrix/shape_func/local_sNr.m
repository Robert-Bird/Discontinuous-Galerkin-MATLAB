function [Nr]=local_sNr(xsi,eta,poly_order)
% Generates shape functions for triangle element K.
% The method comes from: 
% P. Solin, K. Segeth, and I. Dolezel. Higher-order finite element methods. CRC Press, 2003.
%
% Input(s):
% xsi          - Gauss point positions in local coordinate system xsi
% eta          - Gauss point positions in local coordinate system eta
% poly_order   - polynomial order of shape element
%
% Ouput(s):
% Nr - Shape functions

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

% clearing the workspace
% clear;
% clf;

% Order
% Affine coordinates
N1=(eta+1)./2;
N2=-(xsi+eta)./2;
N3=(xsi+1)./2;

% Vertex functions
v1 = N2;
v2 = N3;
v3 = N1;

% Edge functions
Nr=[v1,v2,v3];
c=1;
if poly_order>=2
    N3_N2=kernal_func(N3-N2,poly_order);
    N2_N1=kernal_func(N2-N1,poly_order);
    N1_N3=kernal_func(N1-N3,poly_order);
end
if poly_order>=2
    edge=zeros(max(size(xsi)),3*(poly_order-1));
    for pe=2:poly_order
        edge(:,c)   = N2.*N3.*N3_N2(:,pe-1);
        edge(:,c+1) = N3.*N1.*N1_N3(:,pe-1);
        edge(:,c+2) = N1.*N2.*N2_N1(:,pe-1);
        c=c+3;
    end
    Nr=[v1,v2,v3,edge];
end
% Bubble functions

c=1;
if poly_order>=3
    B=zeros(max(size(xsi)),(poly_order-1)*(poly_order-2)/2);
    for p=3:poly_order
        for n1 = 1:p
            for n2 = 1:p
                if n1+n2<=(p-1) && n1+n2>=(p-1)
                    B(:,c)=(N1.*N2.*N3).*N3_N2(:,n1).*N2_N1(:,n2);
                    c=c+1;
                end
            end
        end
    end
    Nr=[v1,v2,v3,edge,B];
end

end

function [q]=kernal_func(x2,o)
p=o-2;
v = zeros(size(x2,1),p+1);
assert(p>=0,'p must be non-negative');
l = DPolLegrende(x2,p+1);
for i=2:p+2
    v(:,i-1) = -4*sqrt(i-0.5)/(i*(i-1))*l(:,i);
end
q=v;
end
