%RK4 with integrating factor
%set up grid and soliton initial data
clear all;
clc;
N=2^8; dt=.4/N^2; x=(2*pi/N)*(-N/2:N/2-1)';
A=25; B=16; C=36; clf,drawnow
u=3*A^2*sech(.5*(A*(x+2))).^2+3*B^2*sech(.5*(B*(x+1))).^2; %KdV with superposition of solitons
%u=3*B^2*sech(.5*(B*(x+1))).^2; %single soliton %uncomment when necessary

%fourier transform of initial conditions
w=fft(u); 
%k-vector
k=[0:N/2-1 0 -N/2+1:-1]'; 
%linear opeator
ik3=1i*k.^3;

%solve PDE and plot
tmax=0.006; 
nplt=floor((tmax/25)/dt);
nmax=round(tmax/dt);
udata=u;tdata=0;
tic

for n=1:nmax
    t=n*dt;g=-.5i*dt*k;%kdv
    E=exp(dt*ik3/2); E2=E.^2;
    a=g.*fft(real(ifft(w)).^2);
    b=g.*fft(real(ifft(E.*(w+a/2))).^2);
    c=g.*fft(real(ifft(E.*w+b/2)).^2);
    d=g.*fft(real(ifft(E2.*w+E.*c)).^2);
    w=E2.*w+(E2.*a+2*E.*(b+c)+d)/6;
     if mod(n,nplt)==0
          u=real(ifft(w));
          udata=[udata u]; tdata=[tdata t];
     end
end

toc %print elapsed time
nmax %print # of steps involved
waterfall(x,tdata,udata'), view(-20, 25),rotate3d on
xlabel x, ylabel t, axis([-pi pi 0 tmax 0 2000]), grid off
title('Solution to KdV with 1-soliton initial data')