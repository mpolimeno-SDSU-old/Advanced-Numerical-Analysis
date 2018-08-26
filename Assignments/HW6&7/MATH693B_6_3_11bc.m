clear all
tic
h = 1/80;
b = 1;
% mu = 1/h;
mu=5;
tmax = .5;


k = mu*h^2;
td = 0:k:tmax;
q = length(td);
xd = -1:h:1;
p = length(xd);
L = 0:500;

u = zeros(p,1);
for i = 1:p %initial conditions
    u(i,1) = u0(xd(i));
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


alpha = -(b/2)*mu; 
beta = b+b*mu;
mbeta = b-b*mu;

T = zeros(p,p);


    
T(1,1) = beta; %creates T matrix
T(1,2) = 0; 
for i=2:p-1
    T(i,i-1) = alpha;
    T(i,i) = beta;
    T(i,i+1) = alpha;
end
T(p,p) = beta;
T(p,p-1) = 0;

B = zeros(p,p);

B(1,1) = mbeta; %creates B matrix 
B(1,2) = -alpha; 
for i=2:p-1
    B(i,i-1) = -alpha;
    B(i,i) = mbeta;
    B(i,i+1) = -alpha;
end
B(p,p) = mbeta;
B(p,p-1) = -alpha;
    

v = zeros(p,1);
v(1) = w(1,1); %sets initial boundary conditions
for i = 2:p-1
    v(i) = u(i,1); 
end
v(p) = w(1,p);

for i = 2:q
    T(1,1) = ((1-mu)*v(1)+mu/2*v(2))/w(i,1); %changes boundary values
    T(p,p) = (mu/2*v(p-1)+(1-mu)*v(p))/w(i,p);
    v = T\(B*v); %calculates solution vector 
    plot(xd,v,'b-o',xd,w(i,:),'k');
    ylim([0,1.1])
    xlim([-1,1])
    title(['Heat Equation: h = ' num2str(h) ',\mu = ' num2str(mu) ', Time = ' num2str(td(i))])
    xlabel('x')
    ylabel('u')
    grid on
    M(i) = getframe;
end

supnorm = max(abs(w(q,:)-v));

L2norm = norm(w(q,:)-v)

time = toc

