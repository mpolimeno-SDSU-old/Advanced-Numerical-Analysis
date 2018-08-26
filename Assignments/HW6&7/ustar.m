function J = ustar(t,x,p,q) %elle emme
    J = 3/8;
    for jj = 0:p
        for ii = 0:q
            J = J + ((-1)^jj/(pi*(2*jj+1))+2/(pi^2*(2*jj+1)^2))*cos(pi*(2*jj+1)*x)*exp(-pi^2*(2*jj+1)^2)*t...
            + (cos(2*pi*(2*ii+1)*x)/(pi^2*(2*ii+1)^2))*exp(-4*pi^2*(2*ii+1)^2)*t;
        end
    end

end