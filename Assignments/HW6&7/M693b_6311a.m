clear all
tic
h = 1/20;
b = 1;
mu = .4;
tmax = .5;


k = mu*h^2;
t = 0:k:tmax;
qt = length(t);
x = -1:h:(1+h);
px = length(x);
R = 0:200;

for m = 1:px %initial conditions
    u(1,m) = u0(x(m));
end

wp = zeros(qt,px);

for ii = 1:qt %calculates exact solution
    for jj = 1:px
        S1 = ((((-1).^R)./(pi*(2*R+1))) + 2./((pi^2)*((2*R+1).^2)))...
            .*cos(pi*(2*R+1)*x(jj)).*(exp(1).^(-(pi^2)*((2*R+1).^2)*t(ii)));
        S2 = (cos(2*pi*(2*R+1)*x(jj))./((pi^2)*((2*R+1).^2)))...
            .*(exp(1).^(-4*(pi^2)*((2*R+1).^2)*t(ii)));
        wp(ii,jj) = (3/8) + sum(S1) + sum(S2);
    end
end


for i = 1:qt-1 
    u(i,1) = wp(i,1);
    for j = 1:px-2
        u(i+1,j+1) = (1-2*b*mu)*u(i,j+1)+b*mu*(u(i,j+2)+u(i,j));
        u(i,px) = u(i,px-2);
    end
    
end

u(qt,1) = wp(qt,1);

time = toc


for i = 1:qt
    plot(x,u(i,:),'b-o',x,wp(i,:),'k');
    ylim([-0.1,1.5])
    xlim([-1,1])
    title(['Heat equation: \mu = .4, h = ' num2str(h) ', Time = ' num2str(t(i))])
    xlabel('x')
    ylabel('u')
    grid on
    M(i) = getframe;
end


supnorm = max(abs(wp(qt,:)-u(qt,:)))

L2norm = norm(wp(qt,:)-u(qt,:))