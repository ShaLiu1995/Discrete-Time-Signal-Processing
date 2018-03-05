function H_z = zTransform(x)
syms z
H_z = 0;
n = length(x);
for i = 1 : n
    H_z = H_z + x(i) * z^(1-i);
end
display(H_z);