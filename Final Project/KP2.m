%KP2
clear all;
clc;
%set up grid and initial data
N=2^7;
KT=N*2;
L=10;
dt=0.01;
%x = 32*pi*(1:N)'/N;%trefethen paper
Xmesh=linspace(-L,L,KT);
[Xxmesh,Yymesh]=meshgrid(Xmesh,Xmesh); %here where the problem is; have to make is 4096
%u0 = cos(x/16).*(1+sin(x/16)); %trefethen paper iC
%u0 = exp(-x.^2); %trefethen iC
u0 = sin(Xxmesh)*cos(Yymesh).^3; 
w = fft2(u0);%moving to frequency space
clf,drawnow

Dds = 1i.*pi/L*[0:KT/2 -KT/2+1:-1]';
Dy=kron(Dds,ones(KT));
Dx=kron(ones(KT),Dds);
Dx3=Dx.^3;
Dy2=Dx.^2;
%w = reshape(w',KT^2,1);
%Dxn1=[0:1/Dds(1,:)]';


b=1./(Dds(2:end,1));
Dxn1sb=cat(1,zeros(1),b);
Dxn1=(kron(ones(KT),Dxn1sb));
Dx=Dx;
%dt=1/4;% time step

%k=[0:N/2-1 0 -N/2+1:-1]'; 
Lop=Dx3+(Dy2.*Dxn1);
g=-.5i.*Dx.*dt;%Kp2
E=exp(dt.*Lop./2); E2=E.^2;
%solve PDE and plot
uu = u0; tt = 0;
%tmax = 50; nmax = round(tmax/dt); nplt = floor((tmax/75)/dt); %tmax=50
%tmax = 150; nmax = round(tmax/dt); nplt = floor((tmax/100)/dt); %tmax=150
tmax = 10; nmax = round(tmax/dt); nplt = floor((tmax/15)/dt); %tmax=250

nmax=round(tmax/dt);
%h=waitbar(0,'please wait...');
% for n=1:nmax
%     t=n*dt;
%     g=-.5i.*Dx.*dt;%Kp2
%     E=exp(dt.*Lop./2); E2=E.^2;
%     a=Dx.*dt*fft2(real(ifft2(w)).^2);
%     b=Dx.*fft2(real(ifft2(E.*(w+a/2))).^2);
%     c=Dx.*fft2(real(ifft2(E.*w+b/2)).^2);
%     d=Dx.*fft2(real(ifft2(E2.*w+E.*c)).^2);
%     w=E2.*w+(E2.*a+2*E.*(b+c)+d)/6;
%      if mod(n,nplt)==0 
%             u0 = real(ifft2(reshape(w,KT,KT)));waitbar(n/nmax)
%             uu = [uu,u0]; %tt = [tt,t];
%      end
%end
% for n=1:nmax
%     t=n*dt;
%     w = rk4exp(w,dt,g,E,E2,Lop,a,b,c,d);
%     if mod(n,nplt)==0 
%          u0 = real(ifft2(reshape(w,KT,KT)));%waitbar(n/nmax)
%          uu = [uu,u0]; tt = [tt,t];
%     end
% end
for n=1:nmax
    t=n*dt;
    w = rk4exp(w,dt,g,E,E2,Lop,a,b,c,d);
    u0 = real(ifft2(reshape(w',KT^2,1)));%waitbar(n/nmax)
    uu = [uu,u0]; tt = [tt,t];
end

surf(Xxmesh,Yymesh,uu), shading interp, lighting phong, axis tight
% view([-90 90]), colormap(winter); set(gca,'zlim',[-5 50])
%light('color',[0 1 1],'position',[1,1,1])%yellow rgb triple [1 1 0]%blue [0 0 1] %cyan[0 1 1]
%close(h)%close waitbar when plotting is done
%material shiny %shiny is good %metal looks like diarrhea
%xlabel t, ylabel x

% function J = rk4exp(w,dt,g,E,E2,Lop,a,b,c,d)
% t=n*dt;
% g=-.5i.*Dx.*dt;%Kp2
% E=exp(dt.*Lop./2); E2=E.^2;
% a=g.*dt*fft2(real(ifft2(w)).^2);
% b=g.*fft2(real(ifft2(E.*(w+a/2))).^2);
% c=g.*fft2(real(ifft2(E.*w+b/2)).^2);
% d=g.*fft2(real(ifft2(E2.*w+E.*c)).^2);
% J=E2.*w+(E2.*a+2*E.*(b+c)+d)/6;
% end
