clear all


lambda = .9;
h = 1/50;
xd = -4:h:1;
k = lambda*h;
p = length(xd);
td = 0:k:2.7;
q = length(td);


for i = 1:p 
    u(1,i) = f_sol(xd(i));
end


for i = 1:q-1 %lax-friendrichs scheme
    u(i,1) = 1;
    for j = 1:p-2
        u(i+1,j+1) = -lambda*((u(i,j+2) - u(i,j))/2) + ((u(i,j+2) + u(i,j))/2);
    end
end


for i = 1:q-1 %run leap frog
    
    for j = 2:p-2
        u(i+2,j+1) = -lambda*(u(i+1,j+2) - u(i+1,j)) + u(i,j+1);
        u(i+1,p) = 0;
    end
end


for i = 1:q
    v(i,1) = 0;
end


for i = 1:q
    for j = 1:p
        v(i,j) = F_341((xd(j) -  td(i)));
    end
end

for i = 1:q
    plot(xd,u(i,:),'b-o',xd,v(i,:),'k-*');
    ylim([0,2])
    xlim([0,1])
    title(['Wave Function: a = 1, h = ' num2str(h) ',\lambda = ' num2str(lambda) ', Time = ' num2str(td(i))])
    xlabel('x')
    ylabel('u')
    grid on
    M(i) = getframe;
end


for i = 1:q
    e(i,:) = abs((v(i,:)-u(i,:)));
    err(i) = sum(e(i,:)).^2;
end

figure()
plot(td,err,'b-')
title(['SSE vs time: a = 1, h = ' num2str(h) ',\lambda = ' num2str(lambda)])
xlabel('Time')
ylabel('SSE')
grid on




