%one-way wave equation
function J = hwe(x)

if abs(x) <= 1
    J = 1-abs(x);
else 
    J = 0;
end


end