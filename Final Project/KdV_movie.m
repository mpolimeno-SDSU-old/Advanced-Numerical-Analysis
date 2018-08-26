%RK4 with integrating factor
%set up grid and 2-soliton initial data
clear all;
clc;
N=2^8; dt=.4/N^2; x=(2*pi/N)*(-N/2:N/2-1)';
A=25; B=16;clf,drawnow
u=3*A^2*sech(.5*(A*(x+2))).^2+3*B^2*sech(.5*(B*(x+1))).^2; %KdV with superposition of solitons
%u=cos(pi*x);
w=fft(u); k=[0:N/2-1 0 -N/2+1:-1]'; 
ik3=1i*k.^3;
%solve PDE and plot
tmax=0.004; %kdv and burgers
nplt=floor((tmax/25)/dt);
nmax=round(tmax/dt);
%udata=u;tdata=0;%h=waitbar(0,'please wait...');

v = VideoWriter('Kdv_moviewaterfall.avi');
open(v)
for n = 1:nmax
    t=n*dt;g=-.5i*dt*k;%kdv
    E=exp(dt*ik3/2); E2=E.^2;
    a=g.*fft(real(ifft(w)).^2);
    b=g.*fft(real(ifft(E.*(w+a/2))).^2);
    c=g.*fft(real(ifft(E.*w+b/2)).^2);
    d=g.*fft(real(ifft(E2.*w+E.*c)).^2);
    w=E2.*w+(E2.*a+2*E.*(b+c)+d)/6;
    u=real(ifft(w));
    plot(x,u,'linewidth',2.5)
    axis([-pi pi 0 2000])
    M(n) = getframe;
    %%this movie is cooler but was not able to save the file
    %udata=[udata u]; tdata=[tdata t]; 
    %pcolor(x,tdata,udata')
    %shading flat; colormap('jet');
    %title(['KdV-Solitons superposition, t=' num2str(tdata(n))])
    %writeVideo(v,M(n));
end
close(v)