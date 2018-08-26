%Kuramoto-Sivashinky equation
%it evolves into chaos
clear all;
clc;
%set up grid and initial data
N = 2^8;
x=(32*pi/N)*(-N/2:N/2-1)';%periodic grid
%u0 = cos(x/16).*(1+sin(x/16)); %paper IC
u0 = exp(-x.^2); %trefethen IC
w = fft(u0);
clf,drawnow

dt=1/40;
k = [0:N/2-1 0 -N/2+1:-1]'/16; 
k2 = k.^2;%kuramoto
k4 = k.^4; %kuramoto
Lop = (k2-k4);%kuramoto LOP obtained via FFT

%solve PDE and plot
uu = u0; tt = 0;
tmax = 250; nmax = round(tmax/dt); nplt = floor((tmax/375)/dt);

nmax=round(tmax/dt);
h=waitbar(0,'please wait...');
for n=1:nmax
    t=n*dt;
    g=-.5i*k*dt;%kuramoto
    E=exp(dt*Lop/2); E2=E.^2;
    a=g.*fft(real(ifft(w)).^2);
    b=g.*fft(real(ifft(E.*(w+a/2))).^2);
    c=g.*fft(real(ifft(E.*w+b/2)).^2);
    d=g.*fft(real(ifft(E2.*w+E.*c)).^2);
    w=E2.*w+(E2.*a+2*E.*(b+c)+d)/6;
     if mod(n,nplt)==0 
          u0 = real(ifft(w));
          uu = [uu,u0]; tt = [tt,t];
     end
end

surf(tt,x,uu), shading interp, lighting phong, axis tight
view([-90 90]), colormap(winter); set(gca,'zlim',[-5 50])
light('color',[0 1 1],'position',[1,1,1]) %cyan[0 1 1]
material shiny
xlabel t, ylabel x
title('Kuramoto-Sivashinky Equation for $u_{0}=e^{-x^{2}}$','interpreter','latex');