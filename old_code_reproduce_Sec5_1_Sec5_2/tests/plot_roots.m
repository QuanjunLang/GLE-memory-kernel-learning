function plot_roots(r)

figure;
scatter(real(r), imag(r), 130, '.'); hold on;circle(0, 0, 1)
xrange = 1.4;
xlim([-xrange, xrange]);ylim([-xrange, xrange])
title('Prony modes')

end

