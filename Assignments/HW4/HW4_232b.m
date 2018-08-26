%Forward-time central-space scheme 2.3.2b
clear all
clc

lambda = .8;
h = 1/10;
x = -1:h:3;
p = length(x);
k = lambda*h;
t = 0:k:4;
q = length(t);


[X,Y] = meshgrid(x,t);

u = zeros(1,p);
for i = 1:p
    if abs(x(i)) <= 1
        u(1,i) = 1 - abs(x(i));
    else 
        u(1,i) = 0;
    end
end


%left boundary for part b
for i = 1:q
    u(i,1) = -sin(1+i);
end

%run the scheme
for i = 1:q-1
    for j = 1:p-2
        u(i+1,j+1) = -lambda*((u(i,j+2) - u(i,j)))/2 + u(i,j+1);
    end
end

for i = 1:q-1
    for j = 1:p-1
         v(i+1,j+1) = hweb((x(j) -  t(i))); %hweb is for sin(x)
    end
end

u(:,p) = u(:,p-1);

for i = 1:q
    plot(x,u(i,:),'b-o',x,v(i,:),'k-*') 
    ylim([-1,1])
    xlim([-1,3])
    title('FTCS for $u_{0}=\sin(x)$', 'Interpreter', 'latex')
    M(i) = getframe;
end

for i = 1:q
    E(i,:) = abs((v(i,:)-u(i,:)));
    err(i) = sum(E(i,:)).^2;
end
disp(err(i));

figure()
plot(t,err,'b-')
title(['SSE vs time: h = ' num2str(h) ',\lambda = ' num2str(lambda)])
xlabel('Time')
ylabel('SSE')
grid on