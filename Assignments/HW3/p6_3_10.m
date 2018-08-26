clear all
clc

L = 0:500;

h = 1/40;
b = 1;
mu = 10;
tmax = .5;
td = 0:h:tmax;
q = length(td);


k = mu*h^2;
xd = -1:k:1;
p = length(xd);

w = zeros(1,q);
for j = 1:q
    for k = 1:p
        S = ((-1).^L).*((cos(pi*(2*L+1)*xd(k)))./(pi*(2*L+1)))...
            .*(exp(1).^(-(pi^(2))*(2*L+1).^2*td(j)));
        w(j,k) = .5 + 2*sum(S);
    end
end

for i = 1:q
    plot(xd,w(i,:),'k-');
    ylim([0,1])
    xlim([-1,1])
    xlabel('x')
    ylabel('u')
    grid on
    M(i) = getframe;
end


            