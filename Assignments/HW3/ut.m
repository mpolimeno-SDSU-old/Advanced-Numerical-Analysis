function J = ut(t,x,elle)
    J = 1/2;
    for j = 0:elle
        J = J + 2*(-1)^j * ((cos(pi*(2*j+1)*x))/(pi*(2*j+1)))*exp((-pi^2)*(2*j+1)^2*t);
    end

end