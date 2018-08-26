%KP2
clear all;
clc;
%set up grid and initial data
N=2^6;%modes
KT=N*2;
L=10;
dt=0.01;
Xmesh=linspace(-L,L,KT);
[Xxmesh,Yymesh]=meshgrid(Xmesh,Xmesh); 
Xxxmesh=pi./(L*Xxmesh);
Yyymesh=pi./(L*Yymesh);
u0 = sin(Xxxmesh)*cos(Yyymesh).^3;
clf,drawnow

%kronecker tensor products
Dds = 1i.*pi/L*[0:KT/2 -KT/2+1:-1]';
Dy=kron(Dds,ones(KT,1)); %ones(KT) is not the same as np.ones(KT). 
Dx=kron(ones(KT,1),Dds); %ones(KT) creates a KTxKT array in matlab. and a vector of length KT in numpy
Dx3=Dx.^3;
Dy2=6.*Dy.^2;

Dds1 = length(Dds(2:end));
b = ones(Dds1,1)./Dds(2:end);
Dxn1sb = [0; b];
Dxn1=kron(ones(KT,1),Dxn1sb);
Dx=(3/2).*Dx;
 
Lop=Dx3+(Dy2.*Dxn1);%linear operator
g=-.5i.*Dx.*dt;%nonlinearity
E=exp(dt.*Lop./2);%integrating factor (exponential)
%solve PDE and plot
tmax = 10; nmax = round(tmax/dt);

nmax=round(tmax/dt);

Xxxmesh=pi./(L*Xxmesh); %variable change
Yyymesh=pi./(L*Yymesh);

u0 = sin(Xxxmesh)*cos(Yyymesh).^3; %IC's
%reshape and fourier transform initial data
w = reshape(fft2(u0)',KT^2,1);

%call runge kutta 
for n=1:nmax
    t=n*dt;
    w = rk4exp2(w,dt,g,E,KT);
    wnp1 = (real(ifft2(reshape(w.',KT,KT).')))';
end 

pcolor(Xxmesh,Yymesh,wnp1), rotate3d on 
colorbar
xlabel x, ylabel y

