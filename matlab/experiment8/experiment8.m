% experiment 8

list_len = 23;

% parameters
mu_list = linspace(0.1,0.23,list_len); % magnitude of signals
p_list = 50; % parameter size
n_list = 1000000; % sample size
r_list = 1500; % rotation dimension
c0 = 0.25; % sparsity coefficient
s0_list = floor(c0.*p_list); % sparsity

eps = 1; % privacy parameter
del = 1e-2; % privacy parameter
b = 1; % bound on entries of X
gamma = 1/3; %variance of X

lambda = 0.03; % l1 penalty
n_iter = 100;% number of iteration for knockoff power and fdr estimation
q = 0.1; % target fdr

sigma = 1; % noise variance

% results
power_list = zeros(1,length(list_len));
fdr_list = zeros(1,length(list_len));
fdr_th_list = zeros(1,length(list_len));
power_th_list = zeros(1,length(list_len));
fdr_est_th_list = zeros(1,length(list_len));

% Experiment
for i = 1:list_len
    T = DPknockoffs_threshold(n_list,p_list,mu_list(i),r_list,eps,del,b,lambda,gamma,s0_list,q);
    [mean_power,fdr,power_th,fdr_th,fdr_est_th] = DPknockoffs_ker(n_list,p_list,s0_list,mu_list(i),r_list,eps,del,b,n_iter,T,lambda,gamma);
    power_list(i) = mean_power;
    fdr_list(i) = fdr;
    fdr_th_list(i) = fdr_th;
    power_th_list(i) = power_th;
    fdr_est_th_list(i) = fdr_est_th;

end

save("DPknockoffs_expm_8.mat")