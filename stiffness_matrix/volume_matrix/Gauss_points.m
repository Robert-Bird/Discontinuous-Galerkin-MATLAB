function [eta_1,eta_2,w]=Gauss_points(poly_order)
% Generates Gauss point locations for volumetrix integral
%
% Input(s):
% poly_order - Polynomial order
%
% Ouput(s):
% eta_1 - xsi locations of Gauss point for a face
% eta_2 - eta locations of Gauss point for a face
% w     - corresponding wieght of Gauss point

%  Copyright (C) 2017 Robert Bird 
%  $Revision: 1.0 $Date: 2017/06/11 17:09:20 $

N=ceil((poly_order+1));                                                  % Number of Gauss points (could make this less)                                                                                    
v=[-1 -1;1 -1;-1 1];
n=1:N;  nnk=2*n+1; A=[1/3 repmat(1,1,N)./(nnk.*(nnk+2))];
n=2:N; nnk=nnk(n); B1=2/9; nk=n+1; nnk2=nnk.*nnk;
B=4*(n.*nk).^2./(nnk2.*nnk2-nnk2); ab=[A' [2; B1; B']]; s=sqrt(ab(2:N,2));
[V,X]=eig(diag(ab(1:N,1),0)+diag(s,-1)+diag(s,1));
[X,I]=sort(diag(X)); x=(X+1)/2; wx=ab(1,2)*V(1,I)'.^2/4;

N=N-1; N1=N+1; N2=N+2;  y=cos((2*(N:-1:0)'+1)*pi/(2*N+2));
L=zeros(N1,N2);  y0=2;  iter=0;
while max(abs(y-y0))>eps    
    L(:,1)=1;    L(:,2)=y;   
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    y0=y;    y=y0-L(:,N2)./Lp;  iter=iter+1;
end

cd=[ 1, 0, 0; -1, 0, 1; 0, 1,-1]*v; 
t=(1+y)/2;  Wx=abs(det(cd(2:3,:)))*wx;  Wy=1./((1-y.^2).*Lp.^2)*(N2/N1)^2;
[tt,xx]=meshgrid(t,x); yy=tt.*xx;
X=cd(1,1)+cd(2,1)*xx+cd(3,1)*yy;    Y=cd(1,2)+cd(2,2)*xx+cd(3,2)*yy;
w=Wx*Wy';
eta_1=reshape(X,size(X,1)*size(X,2),1);
eta_2=reshape(Y,size(Y,1)*size(Y,2),1);
w=reshape(w,size(w,1)*size(w,2),1);


