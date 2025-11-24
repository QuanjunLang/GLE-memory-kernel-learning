function I = system_settings()
% Settings of learning memory kernel in GLE


%% setting of the DE system and its integrator
I.N         = 501;     % number of x grids
I.L         = 10;       % range of x grids
I.knot_num  = 50;
I.deg       = 3;
I.free_bdry = 1;
I.obs_std   = 1e-3;



I.example_type	= 'mlf';    % {mlf, ....}  
I.h_data_type   = 'grid';   % {prony, grid}


I.obs_dx            = 0.2;
I.prony_p           = 20;
I.prony_N           = 40;
I.prony_obs_std     = 1e-4;


I = update_I(I);

%% h and gamma
plotON = 0;
if plotON
%%
figure;
subplot(131);hold on;
plot(I.xgrid, I.h(I.xgrid));
title('VACF');xlabel('t');
subplot(132);
plot(I.xgrid, I.gamma(I.xgrid));
title('Memory kernel \gamma');xlabel('t');
subplot(133);

loglog(I.log_xgrid, I.gamma(I.log_xgrid));
title('Memory kernel, log scale');xlabel('t');
end




%% save directory
% I.SAVE_DIR = [getenv('HOME'),'/IPS_Graph_data/'];
% if ~exist(I.SAVE_DIR,'dir'); mkdir(I.SAVE_DIR); end
% I.SAVE_DIR_fig = [getenv('HOME'),'/IPS_Graph_data/figures/'];
% if ~exist(I.SAVE_DIR_fig,'dir'); mkdir(I.SAVE_DIR_fig); end
%
% %% Generate the folder for saving figures for paper
% mydir  = pwd;
% idcs   = strfind(mydir,'/');
% newdir = mydir(1:idcs(end)-1);  % This is the main folder for the project
% I.PAPER_FIG_DIR = [newdir, '/Notes/figures'];

end



