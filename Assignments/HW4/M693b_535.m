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

%initial data
u = zeros(1,m); 
for i = 1:m
    if abs(x(i)) <= 1
        u(1,i) = cos(csi*x(i))*cos(0.5*pi*x(i)).^2;
    else
        u(1,i) = 0;
    end
end

%scheme without dissipation
for i = 1:n-1 %leap frog
    u(i+1,1) = 0; %left boundary condition
    for j = 2:m-2
        u(i+2,j+1) = -l*(u(i+1,j+2) - u(i+1,j)) + u(i,j+1);
        u(i+1,m) = u(i,m-1); %right boundary condition %quasi-characteristic extrapolation
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
