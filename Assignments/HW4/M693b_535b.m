clear all
clc

%initialize parameters
l = 0.095;
h = 0.05;
x = -1:h:9;
k = l*h;
t = 0:k:7.5;
m = length(x);
n = length(t);
csi = 5*pi;
a = 1;
d = 1;
eps = 0.5;

%initial data
u = zeros(1,m); 
for i = 1:m
    if abs(x(i)) <= 1
        u(1,i) = cos(csi*x(i))*cos(0.5*pi*x(i)).^2;
    else
        u(1,i) = 0;
    end
end

%scheme with dissipation
for i = 1:n-1
    u(i+1,1) = 0; %left boundary
    for j = 2:m-2
        u(i+2,j+1) = u(i,j+1) - 2*k*a*d*u(i+1,j+1) - (eps*(0.5*h*d).^4)*(u(i,j+1));
        u(i+1,m) = u(i,m-1); %right boundary %quasi-characteristic extrapolation
    end
end

for i = 1:n
    v(i,1) = 0;
end

for i = 1:n
    for j = 1:m
         v(i,j) = v_star((x(j) -  t(i)));
    end
end


for i = 1:n
    plot(x,u(i,:),'b-o',x,v(i,:),'r-')
    ylim([-1,1])
    xlim([-1,2])
    legend('Scheme','v^{\ast}')
    M(i) = getframe;
end
