function [result_h, result_g, result_f] = estimate_h_g_prony(I_h, I_f, F_type)
result_h = prony_method(I_h);

switch F_type
    case '0'
        result_g.g = @(x) -result_h.dh(x);
        result_g.dg = @(x) -result_h.ddh(x);
        result_f = @(x) 0*x;
    otherwise
        result_f = prony_method(I_f);
        result_g.g = @(x) -result_h.dh(x) + result_f.h(x);
        result_g.dg = @(x) -result_h.ddh(x) + result_f.dh(x);
%         
        %         result_g.g = @(x) -result_h.dh(x);
        % result_g.dg = @(x) -result_h.ddh(x);
        % result_f = @(x) 0*x;
end

% obs_grid = 0:I_g.obs_dx:(length(I_g.obs_h_grid)-1)*I_g.obs_dx;
%%
% figure;hold on;
% 
% 
% plot(obs_grid, result_g.h(obs_grid));hold on;
% plot(obs_grid, I_g.obs_h_grid)

%%
% lam = result_h.lam;
% r = result_h.r;
% h = result_h.h;
% w = result_h.w;
% p = length(w);
% 
% 
% Prony_g = @(x) 0*x;
% for t = 1:p
%     Prony_g = @(x) Prony_g(x) + lam(t).*w(t).*exp(lam(t)*x);
% end
% Prony_g = @(x) real(Prony_g(x));
% 
% Prony_dg = @(x) 0*x;
% for t = 1:p
%     Prony_dg = @(x) Prony_dg(x) + lam(t).*lam(t).*w(t).*exp(lam(t)*x);
% end
% Prony_dg = @(x) real(Prony_dg(x));
% 
% 
% 
% result.g = Prony_g;
% result.dg = Prony_dg;


end