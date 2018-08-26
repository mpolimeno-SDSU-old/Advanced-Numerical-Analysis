clear all


lambda = .5;
h = 1/10; %change value of h when needed
x = -1:h:3;
k = lambda*h;
p = length(x);
t = 0:k:2.4;
q = length(t);


for m = 1:p 
    u(1,m) = G(x(m));
end


for i = 1:2 %FTCS scheme
    u(i,1) = 0;
    for j = 1:p-2
        u(i+1,j+1) = -lambda*((u(i,j+2) - u(i,j)))/2 + u(i,j+1);
    end
end


for i = 1:q-1 %leap frog
    
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
        v(i,j) = atan((tan(G((x(j) -  t(i)))))/(1+t(i)*(tan(G((x(j) -  t(i)))))));
    end
end

for i = 1:q
    plot(x,u(i,:),'b-o',x,v(i,:),'k-*');
    ylim([-0.1,1.5])
    xlim([-1,3])
    title(['Wave Function: a = 1, h = ' num2str(h) ',\lambda = ' num2str(lambda) ', Time = ' num2str(t(i))])
    xlabel('x')
    ylabel('u')
    grid on
    M(i) = getframe;
end

supnorm = max(abs(v(q,:)-u(q,:)))

L2norm = norm(v(q,:)-u(q,:))





