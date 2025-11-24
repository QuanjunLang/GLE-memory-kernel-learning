function h_conv_phi = Volterra_convolution(h_grid, phi_grid, dx)

N = length(h_grid);

conv_temp = zeros(N, 1);
conv_temp(1) = 0;
for n = 2:N
    val = (h_grid(1)*phi_grid(n) + h_grid(n)*phi_grid(1))/2;
    val = val + sum(h_grid(2:n-1) .* phi_grid(n-1:-1:2));
    conv_temp(n) = val*dx;
end
h_conv_phi = conv_temp;

end