clear all
tic
h = 1/20;
b = 1;
mu = .4;
tmax = .5;


k = mu*h^2;
td = 0:k:tmax;
q = length(td);
xd = -1:h:(1+h);
p = length(xd);
L = 0:500;


for m = 1:p %initial conditions
    u(1,m) = u0(xd(m));
end

w = zeros(q,p);

for ii = 1:q %calculates exact solution
    for jj = 1:p
        S1 = ((((-1).^L)./(pi*(2*L+1))) + 2./((pi^2)*((2*L+1).^2)))...
            .*cos(pi*(2*L+1)*xd(jj)).*(exp(1).^(-(pi^2)*((2*L+1).^2)*td(ii)));
        S2 = (cos(2*pi*(2*L+1)*xd(jj))./((pi^2)*((2*L+1).^2)))...
            .*(exp(1).^(-4*(pi^2)*((2*L+1).^2)*td(ii)));
        w(ii,jj) = (3/8) + sum(S1) + sum(S2);
    end
end


for i = 1:q-1 
    u(i,1) = w(i,1);
    for j = 1:p-2
        u(i+1,j+1) = (1-2*b*mu)*u(i,j+1)+b*mu*(u(i,j+2)+u(i,j));
        u(i,p) = u(i,p-2);
    end
    
end

u(q,1) = w(q,1);

time = toc


for i = 1:q
    plot(xd,u(i,:),'b-o',xd,w(i,:),'k');
    ylim([-0.1,1.5])
    xlim([-1,1])
    title(['Heat equation: \mu = .4, h = ' num2str(h) ', Time = ' num2str(td(i))])
    xlabel('x')
    ylabel('u')
    grid on
    M(i) = getframe;
end


supnorm = max(abs(w(q,:)-u(q,:)))

L2norm = norm(w(q,:)-u(q,:))





