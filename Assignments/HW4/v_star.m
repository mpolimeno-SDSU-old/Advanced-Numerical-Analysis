function J = v_star(x) %solution

csi = 5*pi;
if abs(x) <= 1
    J = cos(csi*x)*cos(0.5*pi*x).^2;
else
    J = 0;
end
end