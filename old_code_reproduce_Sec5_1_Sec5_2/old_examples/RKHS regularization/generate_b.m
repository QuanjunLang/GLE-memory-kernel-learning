function I = generate_b(I)

dx = I.dx;
basis = I.basis_mat{1};
basis_derivative = I.basis_mat{2};
H = I.H;
n = I.n;


b = zeros(n, 1);
for i = 1:n
    b(i) = H(1)*basis(1, i) + basis_derivative(:, i)'*H*dx;
end

I.b = b;
end


