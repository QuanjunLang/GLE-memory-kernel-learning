function plot_h_conv_phi(corr_func, tgrid, ttl)





figure;hold on;grid on
plot(tgrid, corr_func.h, 'DisplayName','h', 'LineWidth',2);
plot(tgrid, -corr_func.dh, 'DisplayName','-dh', 'LineWidth',2);
% plot(tgrid, h_true(tgrid)*h_tgrid(1), 'DisplayName','h true', 'LineWidth',2);
plot(tgrid, corr_func.h_conv_phi, 'DisplayName','\gamma conv h', 'LineWidth',2);
% plot(tgrid, conv_result_true, 'DisplayName','\gamma conv h', 'LineWidth',2);
plot(tgrid, corr_func.g, '--', 'DisplayName','g', 'LineWidth',2);
plot(tgrid, corr_func.f, '--', 'DisplayName','f', 'LineWidth',2);
% plot(tgrid, g_new, '--', 'DisplayName','new g', 'LineWidth',2);
xlim([0, 30])
legend()
title(ttl)


end