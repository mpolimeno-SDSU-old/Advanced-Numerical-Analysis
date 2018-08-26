clear all


lambda = .5;
h = 1/20;
xd = -1:h:3;
k = lambda*h;
p = length(xd);
td = 0:k:2.4;
q = length(td);


for m = 1:p 
    u(1,m) = F(xd(m));
end


for i = 1:2 %forward-time central-space scheme
    u(i,1) = 0;
    for j = 1:p-2
        u(i+1,j+1) = -lambda*((u(i,j+2) - u(i,j)))/2 + u(i,j+1);
    end
end


for i = 1:q-1 %run leap frog
    
    for j = 2:p-2
        u(i+2,j+1) = -lambda*(u(i+1,j+2) - u(i+1,j)) + u(i,j+1) - 2*k*(sin(u(i,j+1)));
        u(i+1,p) = u(i,p-1);
    end
end

for i = 1:q
    v(i,1) = 0;
end


for i = 1:q
    for j = 1:p
        v(i,j) = atan((tan(F((xd(j) -  td(i)))))/(1+td(i)*(tan(F((xd(j) -  td(i)))))));
    end
end

for i = 1:q
    plot(xd,u(i,:),'b-o',xd,v(i,:),'k-*');
    ylim([-0.1,1.5])
    xlim([-1,3])
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





