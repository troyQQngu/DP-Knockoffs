% experiment 12

list_len1 = 3;
list_len2 = 30;

% parameters
mu_list = 0.13; % magnitude of signals
p_list = 80; % parameter size
n_list = floor(linspace(50000,500000,list_len2)); % sample size
r_list = 1000; % rotation dimension
c0 = 0.25; % sparsity coefficient
s0_list = floor(c0.*p_list); % sparsity

eps_list = linspace(0.2,1,list_len1); % privacy parameter
del = 1e-2; % privacy parameter
b = 1; % bound on entries of X
gamma = 1/3; %variance of X

lambda = 0.03; % l1 penalty
n_iter = 100;% number of iteration for knockoff power and fdr estimation
T_list = 0.05;
sigma = 1; % noise variance

% results
power_list = zeros(list_len1,list_len2);
fdr_list = zeros(list_len1,list_len2);
fdr_th_list = zeros(list_len1,list_len2);
power_th_list = zeros(list_len1,list_len2);
fdr_est_th_list = zeros(list_len1,list_len2);

% Experiment
for i = 1:list_len1 % iterate eps
    for j = 1:list_len2 % iterate n
        %T = DPknockoffs_threshold(n_list,p_list,mu_list,r_list,eps_list(i),del,b,lambda,gamma,s0_list,0.1);
        [mean_power,fdr,power_th,fdr_th,fdr_est_th] = DPknockoffs_ker(n_list(j),p_list,s0_list,mu_list,r_list,eps_list(i),del,b,n_iter,T_list,lambda,gamma);
        power_list(i,j) = mean_power;
        fdr_list(i,j) = fdr;
        fdr_th_list(i,j) = fdr_th;
        power_th_list(i,j) = power_th;
        fdr_est_th_list(i,j) = fdr_est_th;
    end
end

save("DPknockoffs_expm_12.mat")