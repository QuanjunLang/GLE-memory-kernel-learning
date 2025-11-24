function I = generate_b_true(I)


% using direct computation with true gamma
dx = I.dx;
basis = I.basis_mat{1};
basis_derivative = I.basis_mat{2};
H = I.H;
n = I.n;
K = I.K;
gamma_grid = I.gamma_grid;


b = zeros(n, 1);
for i = 1:n
    b(i) = basis(:, i)'*K*gamma_grid'*dx*dx;
end

I.b_true = b;

end