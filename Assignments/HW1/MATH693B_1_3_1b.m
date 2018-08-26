%Forward-time central-space scheme
clear all
clc

lambda = .8;
h = 1/40;
xd = -1:h:3;
p = length(xd);
k = lambda*h;
td = 0:k:2.4;
q=(length(td));
ud = zeros(length(td),1);


[X,Y] = meshgrid(xd,td);

for i = 1:p
    if abs(xd(i)) <= .5
        u(1,i) = cos(pi*xd(i))^2;
    else 
        u(1,i) = 0;
    end
end

for i = 1:q
    u(i,1) = 0;
end


for i = 1:q-1
    for j = 1:p-2
        u(i+1,j+1) = -lambda*((u(i,j+2) - u(i,j)))/2 + u(i,j+1);
    end
end

for i = 1:q-1
    for j = 1:p-1
         v(i+1,j+1) = hwe((xd(j) -  td(i)));
    end
end

u(:,p) = u(:,p-1);


for i = 1:q-90
    plot(xd,u(i,:),'b-o',xd,v(i,:),'k-*');
    ylim([-0.05,1])
    xlim([-1,2])
    M(i) = getframe;
end

title('Forward-time central-space scheme for $h=\frac{1}{40}$ and $\lambda=0.8$', 'Interpreter', 'latex')

