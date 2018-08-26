clear all
clc

h = 1/40;
b = 1;
mu = 10;
tmax = .5;
td = 0:h:tmax;
q = length(td);


k = mu*h^2;
xd = -1:k:1;
p = length(xd);
L = 0:5001;

u = zeros(1,p);
for i = 1:p
    u(i,1) = u0(xd(i));
end

w = zeros(q,p);
for ii = 1:q
    for jj = 1:p
        S = ((-1).^L).*((cos(pi*(2*L+1)*xd(jj)))./(pi*(2*L+1)))...
            .*(exp(1).^(-(pi^(2))*(2*L+1).^2*td(ii)));
        w(ii,jj) = .5 + 2*sum(S);
    end
end

alpha = -(b/2)*mu;
beta = 1+b*mu;
mbeta = 1-b*mu;

T = zeros(p,p);


    
T(1,1) = beta; 
T(1,2) = alpha; 
for i=2:p-1
    T(i,i-1) = alpha;
    T(i,i) = beta;
    T(i,i+1) = alpha;
end
T(p,p) = beta;
T(p,p-1) = alpha;

w = w';

b = zeros(1,q);
for i= 1:q
     b(1) = mbeta*u(1,i) - alpha*u(2,i);
    for j=2:p-1
        b(j) = -alpha*u(j-1,i) + mbeta*u(j,i) + (-alpha)*u(j+1,i);
    end
     b(p) = mbeta*u(p,i) - alpha*u(p-1,i);
    u(:,i+1) = b*inv(T);
end



for i = 1:q
    plot(xd,u(:,i),'b');
    ylim([0,1.1])
    xlim([-1,1])
    title(['Wave Function: a = 1, h = ' num2str(h) ',\mu = ' num2str(mu) ', Time = ' num2str(td(i))])
    xlabel('x')
    ylabel('u')
    grid on
    M(i) = getframe;
end

