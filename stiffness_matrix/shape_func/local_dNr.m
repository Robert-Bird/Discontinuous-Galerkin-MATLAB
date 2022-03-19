function [dNr]=local_dNr(xsi,eta,poly_order)
% Generates shape function derivatives for triangle element K.
% The method comes from: 
% P. Solin, K. Segeth, and I. Dolezel. Higher-order finite element methods. CRC Press, 2003.
%
% Input(s):
% xsi          - Gauss point positions in local coordinate system xsi
% eta          - Gauss point positions in local coordinate system eta
% poly_order   - polynomial order of shape element
%
% Ouput(s):
% dNr - Shape function derivatives

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

% clearing the workspace
% clear;
% clf;
% xsi=-1:0.1:1;
% eta=0;
% Affine coordinates
A1=(eta+1)./2;
A2=-(xsi+eta)./2;
A3=(xsi+1)./2;

% Affine derivatives with respect to local coordinates
A1_x =    0;  A1_e = 1/2;
A2_x = -1/2;  A2_e = -1/2;
A3_x =  1/2;  A3_e =  0;

% Vertex functions
v1=repmat([A2_x;A2_e],max(size(xsi)),1);
v2=repmat([A3_x;A3_e],max(size(xsi)),1);
v3=repmat([A1_x;A1_e],max(size(xsi)),1);
if poly_order>=2
    % Kernal function matrix setup d=0;
    A3_A2_d0=kernal_func((A3-A2),poly_order,0);
    A1_A3_d0=kernal_func((A1-A3),poly_order,0);
    A2_A1_d0=kernal_func((A2-A1),poly_order,0);
    
    % Kernal function matrix setup d=1;
    A3_A2_d1=kernal_func((A3-A2),poly_order,1);
    A1_A3_d1=kernal_func((A1-A3),poly_order,1);
    A2_A1_d1=kernal_func((A2-A1),poly_order,1);
end

% Edge functions
dNr=[v1,v2,v3];
c=1;
if poly_order>=2
    edge=zeros(max(size(xsi))*2,3*(poly_order-1));
    for pe=2:poly_order
        % differentiate with xsi
        edge(1:2:end,c)   = A2_x.*A3.*A3_A2_d0(:,pe-1) +A2.*A3_x.*A3_A2_d0(:,pe-1) +A2.*A3.*A3_A2_d1(:,pe-1);
        edge(1:2:end,c+1) = A3_x.*A1.*A1_A3_d0(:,pe-1) +A3.*A1_x.*A1_A3_d0(:,pe-1) +A3.*A1.*A1_A3_d1(:,pe-1).*(-1/2);
        edge(1:2:end,c+2) = A1_x.*A2.*A2_A1_d0(:,pe-1) +A1.*A2_x.*A2_A1_d0(:,pe-1) +A1.*A2.*A2_A1_d1(:,pe-1).*(-1/2);
        % differentiate with eta
        edge(2:2:end,c)   = A2_e.*A3.*A3_A2_d0(:,pe-1) +A2.*A3_e.*A3_A2_d0(:,pe-1) +A2.*A3.*A3_A2_d1(:,pe-1).*(1/2);
        edge(2:2:end,c+1) = A3_e.*A1.*A1_A3_d0(:,pe-1) +A3.*A1_e.*A1_A3_d0(:,pe-1) +A3.*A1.*A1_A3_d1(:,pe-1).*(1/2);
        edge(2:2:end,c+2) = A1_e.*A2.*A2_A1_d0(:,pe-1) +A1.*A2_e.*A2_A1_d0(:,pe-1) +A1.*A2.*A2_A1_d1(:,pe-1).*(-1);
        c=c+3;
    end
    dNr=[v1,v2,v3,edge];
end
% Bubble functions
c=1;
if poly_order>=3
    B=zeros(max(size(xsi))*2,(poly_order-1)*(poly_order-2)/2);
    for p=3:poly_order
        for n1 = 1:p
            for n2 = 1:p
                if n1+n2<=(p-1) && n1+n2>=(p-1)
                    B(1:2:end,c)=(A1_x.*A2.*A3).*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...
                        +(A1.*A2_x.*A3).*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...
                        +(A1.*A2.*A3_x).*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...
                        +(A1.*A2.*A3)  .*A3_A2_d1(:,n1).*A2_A1_d0(:,n2)...
                        +(A1.*A2.*A3)  .*A3_A2_d0(:,n1).*A2_A1_d1(:,n2).*(-1/2);
                    
                    B(2:2:end,c)=(A1_e.*A2.*A3).*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...
                        +(A1.*A2_e.*A3).*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...
                        +(A1.*A2.*A3_e).*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...
                        +(A1.*A2.*A3)  .*A3_A2_d1(:,n1).*A2_A1_d0(:,n2).*(1/2)...
                        +(A1.*A2.*A3)  .*A3_A2_d0(:,n1).*A2_A1_d1(:,n2).*(-1);
                    c=c+1;
                end
            end
        end
    end
    dNr=[v1,v2,v3,edge,B];
end

end


% Differentiated Kernal function
function [q]=kernal_func(x2,o,d)
if d==0
    p=o-2;
    v = zeros(size(x2,1),p+1);
    assert(p>=0,'p must be non-negative');
    l = DPolLegrende(x2,p+1);
    for i=2:p+2
        v(:,i-1) = -4*sqrt(i-0.5)/(i*(i-1))*l(:,i);
    end
    q=v;
end

if d==1
    p=o-2;
    v = zeros(size(x2,1),p+1);
    assert(p>=0,'p must be non-negative');
    l = DDPolLegrende(x2,p+1);
    for i=2:p+2
        v(:,i-1) = -4*sqrt(i-0.5)/(i*(i-1))*l(:,i);
    end
    q=v;
end
end
