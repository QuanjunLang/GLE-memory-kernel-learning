function I = generate_A(I)

dx = I.dx;
basis_mat = I.basis_mat;

basis = basis_mat{1};

A = basis'*I.K*basis*dx*dx;

I.A = A;

end


