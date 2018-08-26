function J = G(x)

if abs(x) <= .5
    J = cos(pi*x)^2;
else 
    J = 0;
end
end