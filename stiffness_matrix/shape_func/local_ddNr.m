function [ddNr]=local_ddNr(xsi,eta,pb)
% Generates second order derivatives for triangle element K.
% The method comes from: 
% P. Solin, K. Segeth, and I. Dolezel. Higher-order finite element methods. CRC Press, 2003.
%
% Input(s):
% xsi  - Gauss point positions in local coordinate system xsi
% eta  - Gauss point positions in local coordinate system eta
% pb   - polynomial order of shape element
%
% Ouput(s):
% ddNr - Shape function derivatives

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

% Affine derivatives with respect to local coordinates, all higher order
% derivatives are zero
A1_x =    0;  A1_e = 1/2;
A2_x = -1/2;  A2_e = -1/2;
A3_x =  1/2;  A3_e =  0;

% Vertex functions second order derivatives
v1=repmat([0;0],max(size(xsi)),1);
v2=repmat([0;0],max(size(xsi)),1);
v3=repmat([0;0],max(size(xsi)),1);

if pb>=2
    % Matrices of Kernal function values
    A3_A2_d0=kernal_func((A3-A2),pb,0);
    A3_A2_d1=kernal_func((A3-A2),pb,1);
    A3_A2_d2=kernal_func((A3-A2),pb,2);
    
    % Matrices of Kernal function values
    A1_A3_d0=kernal_func((A1-A3),pb,0);
    A1_A3_d1=kernal_func((A1-A3),pb,1);
    A1_A3_d2=kernal_func((A1-A3),pb,2);
    
    % Matrices of Kernal function values
    A2_A1_d0=kernal_func((A2-A1),pb,0);
    A2_A1_d1=kernal_func((A2-A1),pb,1);
    A2_A1_d2=kernal_func((A2-A1),pb,2);
end
% Edge functions
ddNr=[v1,v2,v3];
c=1;
if pb>=2
    edge=zeros(max(size(xsi))*2,3*(pb-1));
    for pe=2:pb
        % differentiate with xsi(x)
        % Edge 1 - some terms are zero, so only terms additions not 10
        edge(1:2:end,c) = A2_x.*A3_x.*A3_A2_d0(:,pe-1)...         % first line
            +A2_x.*A3  .*A3_A2_d1(:,pe-1)...         %
            +A2_x.*A3_x.*A3_A2_d0(:,pe-1)...         % Second line
            +A2  .*A3_x.*A3_A2_d1(:,pe-1)...         %
            +A2_x.*A3  .*A3_A2_d1(:,pe-1)...         % Third line
            +A2  .*A3_x.*A3_A2_d1(:,pe-1)...         %
            +A2  .*A3  .*A3_A2_d2(:,pe-1);           % Second order kernal function derivative, note: d(A3-A2)/dx=1
        edge(1:2:end,c+1)=A3_x.*A1_x.*A1_A3_d0(:,pe-1)...         % first line
            +A3_x.*A1  .*A1_A3_d1(:,pe-1).*(-1/2)... %
            +A3_x.*A1_x.*A1_A3_d0(:,pe-1)...         % Second line
            +A3  .*A1_x.*A1_A3_d1(:,pe-1).*(-1/2)... %
            +A3_x.*A1  .*A1_A3_d1(:,pe-1).*(-1/2)... % Third line
            +A3  .*A1_x.*A1_A3_d1(:,pe-1).*(-1/2)... %
            +A3  .*A1  .*A1_A3_d2(:,pe-1).*(-1/2)^2; % Second order kernal function derivative
        edge(1:2:end,c+2)=A1_x.*A2_x.*A2_A1_d0(:,pe-1)...         % first line
            +A1_x.*A2  .*A2_A1_d1(:,pe-1).*(-1/2)... %
            +A1_x.*A2_x.*A2_A1_d0(:,pe-1)...         % Second line
            +A1  .*A2_x.*A2_A1_d1(:,pe-1).*(-1/2)... %
            +A1_x.*A2  .*A2_A1_d1(:,pe-1).*(-1/2)... % Third line
            +A1  .*A2_x.*A2_A1_d1(:,pe-1).*(-1/2)... %
            +A1  .*A2  .*A2_A1_d2(:,pe-1).*(-1/2)^2; % Second order kernal function derivative
        % differentiate with eta(e)
        edge(2:2:end,c) =  A2_e.*A3_e.*A3_A2_d0(:,pe-1)...         % first line
            +A2_e.*A3  .*A3_A2_d1(:,pe-1).*(1/2)...  %
            +A2_e.*A3_e.*A3_A2_d0(:,pe-1)...         % Second line
            +A2  .*A3_e.*A3_A2_d1(:,pe-1).*(1/2)...  %
            +A2_e.*A3  .*A3_A2_d1(:,pe-1).*(1/2)...  % Third line
            +A2  .*A3_e.*A3_A2_d1(:,pe-1).*(1/2)...  %
            +A2  .*A3  .*A3_A2_d2(:,pe-1).*(1/2)^2;  % Second order kernal function derivative, note: d(A3-A2)/dx=1
        edge(2:2:end,c+1)=A3_e.*A1_e.*A1_A3_d0(:,pe-1)...         % first line
            +A3_e.*A1  .*A1_A3_d1(:,pe-1).*(1/2)... %
            +A3_e.*A1_e.*A1_A3_d0(:,pe-1)...         % Second line
            +A3  .*A1_e.*A1_A3_d1(:,pe-1).*(1/2)... %
            +A3_e.*A1  .*A1_A3_d1(:,pe-1).*(1/2)... % Third line
            +A3  .*A1_e.*A1_A3_d1(:,pe-1).*(1/2)... %
            +A3  .*A1  .*A1_A3_d2(:,pe-1).*(1/2)^2; % Second order kernal function derivative
        edge(2:2:end,c+2)=A1_e.*A2_e.*A2_A1_d0(:,pe-1)...         % first line
            +A1_e.*A2  .*A2_A1_d1(:,pe-1).*(-1)... %
            +A1_e.*A2_e.*A2_A1_d0(:,pe-1)...         % Second line
            +A1  .*A2_e.*A2_A1_d1(:,pe-1).*(-1)... %
            +A1_e.*A2  .*A2_A1_d1(:,pe-1).*(-1)... % Third line
            +A1  .*A2_e.*A2_A1_d1(:,pe-1).*(-1)... %
            +A1  .*A2  .*A2_A1_d2(:,pe-1).*(-1)^2; % Second order kernal function derivative
        c=c+3;
    end
    ddNr=[v1,v2,v3,edge];
end
% Bubble functions
c=1;
if pb>=3
    B=zeros(max(size(xsi))*2,(pb-1)*(pb-2)/2);
    for p=3:pb
        for n1 = 1:p
            for n2 = 1:p
                if n1+n2<=(p-1) && n1+n2>=(p-1)
                    B(1:2:end,c)=(A1_x.*A2_x.*A3)  .*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...%part 1
                        +(A1_x.*A2  .*A3_x).*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...
                        +(A1_x.*A2  .*A3)  .*A3_A2_d1(:,n1).*A2_A1_d0(:,n2)...
                        +(A1_x.*A2  .*A3)  .*A3_A2_d0(:,n1).*A2_A1_d1(:,n2).*(-1/2)...
                        +(A2_x.*A1_x.*A3)  .*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...%part 2
                        +(A2_x.*A1  .*A3_x).*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...
                        +(A2_x.*A1  .*A3)  .*A3_A2_d1(:,n1).*A2_A1_d0(:,n2)...
                        +(A2_x.*A1  .*A3)  .*A3_A2_d0(:,n1).*A2_A1_d1(:,n2).*(-1/2)...
                        +(A3_x.*A2_x.*A1)  .*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...%part 3
                        +(A3_x.*A2  .*A1_x).*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...
                        +(A3_x.*A2  .*A1)  .*A3_A2_d1(:,n1).*A2_A1_d0(:,n2)...
                        +(A3_x.*A2  .*A1)  .*A3_A2_d0(:,n1).*A2_A1_d1(:,n2).*(-1/2)...
                        ...
                        +(A1_x.*A2.*A3)    .*A3_A2_d1(:,n1).*A2_A1_d0(:,n2)...%part 4
                        +(A1.*A2_x.*A3)    .*A3_A2_d1(:,n1).*A2_A1_d0(:,n2)...
                        +(A1.*A2.*A3_x)    .*A3_A2_d1(:,n1).*A2_A1_d0(:,n2)...
                        +(A3.*A2.*A1)      .*A3_A2_d2(:,n1).*A2_A1_d0(:,n2)... % second order bit
                        +(A1.*A2.*A3)      .*A3_A2_d1(:,n1).*A2_A1_d1(:,n2).*(-1/2)...
                        ...
                        +(A1_x.*A2.*A3)    .*A3_A2_d0(:,n1).*A2_A1_d1(:,n2).*(-1/2)...%part 5
                        +(A1.*A2_x.*A3)    .*A3_A2_d0(:,n1).*A2_A1_d1(:,n2).*(-1/2)...
                        +(A1.*A2.*A3_x)    .*A3_A2_d0(:,n1).*A2_A1_d1(:,n2).*(-1/2)...
                        +(A3.*A2.*A1)      .*A3_A2_d1(:,n1).*A2_A1_d1(:,n2).*(-1/2)...
                        +(A1.*A2.*A3)      .*A3_A2_d0(:,n1).*A2_A1_d2(:,n2).*(-1/2)^2; % second order bit
                    
                    B(2:2:end,c)=(A1_e.*A2_e.*A3)  .*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...%part 1
                        +(A1_e.*A2  .*A3_e).*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...
                        +(A1_e.*A2  .*A3)  .*A3_A2_d1(:,n1).*A2_A1_d0(:,n2).*(1/2)...
                        +(A1_e.*A2  .*A3)  .*A3_A2_d0(:,n1).*A2_A1_d1(:,n2).*(-1)...
                        ...
                        +(A2_e.*A1_e.*A3)  .*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...%part 2
                        +(A2_e.*A1  .*A3_e).*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...
                        +(A2_e.*A1  .*A3)  .*A3_A2_d1(:,n1).*A2_A1_d0(:,n2).*(1/2)...
                        +(A2_e.*A1  .*A3)  .*A3_A2_d0(:,n1).*A2_A1_d1(:,n2).*(-1)...
                        ...
                        +(A3_e.*A2_e.*A1)  .*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...%part 3
                        +(A3_e.*A2  .*A1_e).*A3_A2_d0(:,n1).*A2_A1_d0(:,n2)...
                        +(A3_e.*A2  .*A1)  .*A3_A2_d1(:,n1).*A2_A1_d0(:,n2).*(1/2)...
                        +(A3_e.*A2  .*A1)  .*A3_A2_d0(:,n1).*A2_A1_d1(:,n2).*(-1)...
                        ...
                        +(A1_e.*A2  .*A3)  .*A3_A2_d1(:,n1).*A2_A1_d0(:,n2).*(1/2)...%part 4
                        +(A1  .*A2_e.*A3)  .*A3_A2_d1(:,n1).*A2_A1_d0(:,n2).*(1/2)...
                        +(A1  .*A2  .*A3_e).*A3_A2_d1(:,n1).*A2_A1_d0(:,n2).*(1/2)...
                        +(A3  .*A2  .*A1)  .*A3_A2_d2(:,n1).*A2_A1_d0(:,n2).*(1/2)^2 ... % second order bit
                        +(A1  .*A2  .*A3)  .*A3_A2_d1(:,n1).*A2_A1_d1(:,n2).*(-1/2)...
                        ...
                        +(A1_e.*A2  .*A3)  .*A3_A2_d0(:,n1).*A2_A1_d1(:,n2).*(-1)...%part 5
                        +(A1  .*A2_e.*A3)  .*A3_A2_d0(:,n1).*A2_A1_d1(:,n2).*(-1)...
                        +(A1  .*A2  .*A3_e).*A3_A2_d0(:,n1).*A2_A1_d1(:,n2).*(-1)...
                        +(A3  .*A2  .*A1)  .*A3_A2_d1(:,n1).*A2_A1_d1(:,n2).*(-1/2)...
                        +(A1  .*A2  .*A3)  .*A3_A2_d0(:,n1).*A2_A1_d2(:,n2).*(-1)^2; % second order bit
                    c=c+1;
                end
            end
        end
    end
    ddNr=[v1,v2,v3,edge,B];
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

if d==2
    p=o-2;
    v = zeros(size(x2,1),p+1);
    assert(p>=0,'p must be non-negative');
    l = DDDPolLegrende(x2,p+1);
    for i=2:p+2
        v(:,i-1) = -4*sqrt(i-0.5)/(i*(i-1))*l(:,i);
    end
    q=v;
end
end
