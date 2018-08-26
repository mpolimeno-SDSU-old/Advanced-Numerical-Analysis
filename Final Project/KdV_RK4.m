%set up grid and 2-soliton initial data %RK4 without integrating factor
clear all;
clc;
N=2^8; dt=.4/N^2; x=(2*pi/N)*(-N/2:N/2-1)';
A=25; B=16; C=9; clf,drawnow
u=3*A^2*sech(.5*(A*(x+2))).^2+3*B^2*sech(.5*(B*(x+1))).^2; %KdV %2 solitons
%u=3*A^2*sech(.5*(A*(x+2))).^2; %1 soliton
w=fft(u); 
k=[0:N/2-1 0 -N/2+1:-1]'; 
g=-.5i*dt.*k;
%solve PDE and plot
tmax=0.006; %kdv and burgers;
nplt=floor((tmax/25)/dt);
nmax=round(tmax/dt);
udata=u;tdata=0;%h=waitbar(0,'please wait...');
tic
% for n=1:nmax
%     t=n*dt;
%     w=rk4(w,dt);
%     if mod(n,nplt)==0
%          u=real(ifft(w));waitbar(n/nmax)
%          udata=[udata u]; tdata=[tdata t];
%     end
% end
% waterfall(x,tdata,udata'), view(-20, 25),rotate3d on
% xlabel x, ylabel t, axis([-pi pi 0 tmax 0 2000]), grid off
% set(gca, 'ztick', [0 2000]), close(h), pbaspect([1 1 .13])
for i = 1:nmax
    t=i*dt;
    w=rk4(w,dt);
    u=real(ifft(w));
    udata=[udata u]; tdata=[tdata t];
    pcolor(x,tdata,udata')
    colorbar;
    shading flat; colormap('jet');
    %title(['Advection-Diffusion, t=' num2str(t(i))])
    M(i) = getframe;
%     writeVideo(v,M(i));
end
toc

function J = rk4(w,dt)
a=dt*fft(real(ifft(w)));
b=dt*dt*fft(real(ifft((w+a/2))));
c=dt*dt*fft(real(ifft(w+b/2)));
d=dt*dt*fft(real(ifft(w+c)));
J=w+(dt*(a+2*(b+c)+d))/6;
end
