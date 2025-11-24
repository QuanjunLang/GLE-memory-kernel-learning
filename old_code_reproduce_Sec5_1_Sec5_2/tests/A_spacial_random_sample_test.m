clc
close all
clear all

random_seed = 15;
rng(random_seed)



memory_para.p = 1;
memory_para.q = 3;
memory_para.plotON = 1;
[R, S, R_lap, u_seq, eta_seq, h_true, h_true_lap] = generate_memory_kernel(memory_para);


[v, dt, tgrid, T0, correlated_noise] = generate_trajectories(R, S, F, M, traj_para, random_seed);