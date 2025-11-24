function I = generate_K(I)
dx = I.dx;
N = I.N;
h_grid = I.h_grid;
xgrid = I.xgrid;

H = zeros(N, 1);
h_pad = [h_grid, zeros(size(xgrid))];
for t = 1:N
    H(t) = sum(h_grid.*h_pad(t:t+I.N-1))*dx*dx;
end

K = zeros(N, N);
for i = 1:N
    for j = 1:N
        K(i, j) = H(abs(i-j)+1);
    end
end



I.K = K;
I.H = H;


% figure;
% [XX, YY] = meshgrid(xgrid, xgrid);
% surf(XX, YY, K);
% shading flat
end

