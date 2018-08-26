% Solve the Poisson Equation using SOR method
clear; close all;
L = 1;
dx = .1;
dy = dx;
h= dx;
n = round(L/dx);
m = n;
x = linspace(0,L,n)';
y = linspace(0,L,m)';
T_0 = zeros(n,m);

bb = 2/(1 + pi*h); 

% Boundary Condition cosxsiny
% At x= 0; 
for i =1:n
    T_0(i,1) = cos(x(1))*sin(y(i));  
% At x = 1
    T_0(i,m) = cos(x(m)).*sin(y(i));
end
% At y = 0;
for j =1:m
    T_0(1,j) = cos(x(j))*sin(y(1));
% At y = 1;
    T_0(m,j) = cos(x(j))*sin(y(m));
end

T= T_0;

S = zeros(n,m); 

for i = 1:n
   for j=1:m
         S(i,j)= -2*cos(x(i))*sin(y(j));
    end
end


max_step = 1000;


for l=1:max_step
 for i=2:n-1
     % Boundary Conditions
     
     for j=2:m-1
        T(i,j)=bb*0.25*(T(i+1,j)+...
        T(i,j+1)+T(i-1,j)+T(i,j-1)-h^2*S(i,j))+(1.0-bb)*T(i,j);
        % Boundary Conditions
       
     end
 end
 % find residual
 res=0;
 for i=2:n-1
     for j=2:m-1
        res=res+abs(T(i+1,j)+...
        T(i,j+1)+T(i-1,j)+T(i,j-1)-4*T(i,j))/h^2 - S(i,j);
     end     
     
 end
 if rem(l,2) == 0
     figure(1);
     surf(x,y,T);
     title([' SOR for Poisson Equation, at iteration ',num2str(l),' with'...
        ' h =',num2str(h)]);
     xlabel('x')
     ylabel('y')
     zlabel('z')
     axis([ 0 1 0 1 0 1]);
     getframe(gcf);
 end
 
 l,res/((m-2)*(n-2)) % Print iteration and residual
 if (res/((m-2)*(n-2)) < 0.001)
     break
 end
 disp('Norm:')
 Lnorm = (norm(T)^2*h^2)
end