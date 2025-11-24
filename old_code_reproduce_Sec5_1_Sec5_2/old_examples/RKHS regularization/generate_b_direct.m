function I = generate_b_direct(I)

dx = I.dx;
basis = I.basis_mat{1};
basis_derivative = I.basis_mat{2};
H = I.H;
n = I.n;
N = I.N;
h_grid = I.h_grid;
g_grid = I.g_grid;


g_grid_pad = [g_grid, zeros(1, N)];
HG = zeros(N, 1);
for i = 1:N
    HG(i) = sum(h_grid.*g_grid_pad(i:N+i-1))*dx*dx;
end

b = -HG'*basis*dx;


I.b_direct = b';
end


