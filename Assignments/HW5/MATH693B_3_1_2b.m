clear all

lambda = .8;
h = 1/10; %change the value of h as needed
x = -1:h:1;
k = lambda*h;
p = length(x);
t = 0:k:1.2;
q = length(t);


for k = 1:p 
    u(1,k) = G8(x(k)); %set initial data
end


for i = 1:q-1 %lax-wendroff scheme
    for j = 1:p-2
        u(i,1) = G8((x(j) -  t(i))); %set boundary data1
        u(i+1,j+1) = u(i,j+1) - (lambda/2)*(u(i,j+2)-u(i,j)) + ((lambda^2)/2)*(u(i,j+2) - 2*u(i,j+1) + u(i,j));
    end
end
u(q,1) = G8((x(1) -  t(q)));
u(:,p) = u(:,1);

for i = 1:p %boundary data for exact solution
    v(i,1) = 0;
end


for i = 1:q %calculates exact solution from initial data
    for j = 1:p
        v(i,j) = G8((x(j) -  t(i)));
    end
end

for i = 1:q
    plot(x,u(i,:),'b-o',x,v(i,:),'k-*');
    ylim([-1,1])
    xlim([-1,1])
    title(['Wave Function: a = 1, h = ' num2str(h) ',\lambda = ' num2str(lambda) ', Time = ' num2str(t(i))])
    xlabel('x')
    ylabel('u')
    grid on
    M(i) = getframe;
end


for i = 1:q %calculates Sum of Squared Errors
    e(i,:) = abs((v(i,:)-u(i,:)));
    err(i) = sum(e(i,:)).^2;
end

supnorm = max(abs(v(q,:)-u(q,:)))

L2norm = norm(v(q,:)-u(q,:))

figure()
plot(td,err,'b-') %plots SSEs at all calculated times
title(['SSE vs time: a = 1, h = ' num2str(h) ',\lambda = ' num2str(lambda)])
xlabel('Time')
ylabel('SSE')
grid on
