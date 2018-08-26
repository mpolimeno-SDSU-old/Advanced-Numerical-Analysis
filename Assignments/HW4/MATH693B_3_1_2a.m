clear all

lambda = .8;
h = 1/80;
xd = -1:h:1;
k = lambda*h;
p = length(xd);
td = 0:k:1.2;
q = length(td);


for i = 1:p 
    u(1,i) = F312(xd(i)); %set initial data
end


for i = 1:q %forward-time backward-space scheme
    for j = 1:p-1
        u(i,1) = u(i,p); %set boundary data
        u(i+1,j+1) = (1-lambda)*u(i,j+1) + lambda*u(i,j);
    end
end

u(:,p) = u(:,1);

for i = 1:p %boundary data for exact solution
    v(i,1) = 0;
end


for i = 1:q %calculates exact solution from initial data
    for j = 1:p
        v(i,j) = F312((xd(j) -  td(i)));
    end
end

for i = 1:q
    plot(xd,u(i,:),'b-o',xd,v(i,:),'k-*');
    ylim([-1,1])
    xlim([-1,1])
    title(['Wave Function: a = 1, h = ' num2str(h) ',\lambda = ' num2str(lambda) ', Time = ' num2str(td(i))])
    xlabel('x')
    ylabel('u')
    grid on
    M(i) = getframe;
end

supnorm = max(abs(v(q,:)-u(q,:)))

L2norm = norm(v(q,:)-u(q,:))


for i = 1:q %calculates Sum of Squared Errors
    e(i,:) = abs((v(i,:)-u(i,:)));
    err(i) = sum(e(i,:)).^2;
end

figure()
plot(td,err,'b-') %plots SSEs at all calculated times
title(['SSE vs time: a = 1, h = ' num2str(h) ',\lambda = ' num2str(lambda)])
xlabel('Time')
ylabel('SSE')
grid on
