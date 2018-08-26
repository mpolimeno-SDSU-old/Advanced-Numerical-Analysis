clc
clear all

u0 = @(x) (abs(x) < 0.5)/(1-abs(x))+(abs(x) == 0.5)/4;
tmax = 1/2;
p = 100;
q = 100;
hvals = [1/10 1/20 1/40 1/80];

for lambda_switch = 1:2
    for h = hvals
        if lambda_switch == 1
            mu = 1/h; %lambda=mu*h and for lambda=1 we have mu=1/h
        else 
            mu = 5;
        end
        
        k = mu*h^2;
        x = (-1:h:1)';
        m = length(x);
        b = ones(m,1);
        v = zeros(m,1);
        time = 0;
        
        coeff_matrx_n1 = [(-mu/2)*b (1+mu)*b (-mu/2)*b];
        coeff_matrx_n  = [(mu/2)*b (1-mu)*b (mu/2)*b];
        
        A = spdiags(coeff_matrx_n1,[-1 0 1], m,m); %make triadiagonal matrix %diagonals are at the -1 (meaning 1 down), main and 1 (meaning 1 up)
        B = spdiags(coeff_matrx_n, [-1 0 1], m,m);
        A(1,2) = 0;
        A(m,m-1) = 0;
        
        v(1) = ustar(0,1,p,q);
        for M = 2:m-1
            v(M) = u0(x(M)); %initial data
        end
        v(m) = ustar(0,m,p,q);
        
        while time < tmax
            time = time+k; %compute time
            
            A(1,1) = ((1-mu)*v(1)+mu/2*v(2))/ustar(time,1,p,q);
            A(m,m) = (mu/2*v(m-1)+(1-mu)*v(m))/ustar(time,m,p,q);
            v = A\(B*v);
        end
        
        uxsol = ustar(time,x,p,q);
        figure;
        hold on
        if lambda_switch == 1
            type = '\lambda=1';
        else
            type = '\mu = 5';
        end
        title(['h = ' num2str(h) ',' type]);
        plot(x,uxsol,'*',x,v,'x-');
        xlabel('x')
        ylabel('u')
        axis([-1 1 0 1])
        hold off;
        if lambda_switch == 1
            disp(['h: ' num2str(h)]);
            disp(['Supremum Norm: ' num2str(max(abs(uxsol-v)))]);
            disp(['L2 Norm: ' num2str(norm(uxsol-v))]);
        end
    end
end
     