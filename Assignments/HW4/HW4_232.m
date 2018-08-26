%Forward-time central-space scheme
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

%left boundary for part a
for i = 1:q  
    u(i,1) = 0;
end

%run the scheme
for i = 1:q-1
    for j = 1:p-2
        u(i+1,j+1) = -lambda*((u(i,j+2) - u(i,j)))/2 + u(i,j+1);
    end
end

for i = 1:q-1                             
    for j = 1:p-1
         v(i+1,j+1) = hwe((x(j) -  t(i))); %hwe is for 1-abs(x)
    end
end

u(:,p) = u(:,p-1);


for i = 1:q                              
    plot(x,u(i,:),'b-o',x,v(i,:),'k-*');
    ylim([-0.6,1])
    xlim([-1,2])
    title('FTCS scheme for $u_{0}=1-|x|$ or $0$' , 'Interpreter', 'latex')
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
