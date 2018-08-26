function J = f_sol(x)
if x <= 0
    J = 1;
else
    J = cos(2*pi*x);
end