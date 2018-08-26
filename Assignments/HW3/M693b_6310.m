clc
clear all

u0 = @(x) (abs(x) < 0.5)+(abs(x) == 0.5)/2;

tmax = 1/2;
elle = 100;
hvals = [1/10 1/20 1/40];

for lambda_switch = 1:2
    for h = hvals
        if lambda_switch == 1
            mu = 1/h; %lambda=mu*h and for lambda=1 we have mu=1/h
        else 
            mu = 10;
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
        
        v(1) = ut(0,1,elle);
        for M = 2:m-1
            v(M) = u0(x(M)); %initial data
        end
        v(m) = ut(0,m,elle);
        
        while time < tmax
            time = time+k; %compute time
            
            A(1,1) = ((1-mu)*v(1)+mu/2*v(2))/ut(time,1,elle);
            A(m,m) = (mu/2*v(m-1)+(1-mu)*v(m))/ut(time,m,elle);
            v = A\(B*v);
        end
        
        uxsol = ut(time,x,elle);
        figure;
        hold on
        if lambda_switch == 1
            type = '\lambda=1';
        else
            type = '\mu = 10';
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
        
            
            
            
    
    