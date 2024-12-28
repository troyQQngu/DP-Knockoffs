% experiment 23. Same experiment as 22, but we will show the JLT has better
% performance when signal magnitude is large

list_len = 13;

% parameters
mu_list = 2; % magnitude of signals
p_list = 50; % parameter size
n_list = floor(linspace(10000,120000,list_len)); % sample size
r_list = 10000; % rotation dimension
c0 = 0.3; % sparsity coefficient
s0_list = floor(c0.*p_list); % sparsity

eps = 0.2; % privacy parameter
del =  0.01; % privacy parameter
b = 1; % bound on entries of X
gamma = 1/3; %variance of X

lambda = 0.025; % l1 penalty
n_iter = 100;% number of iteration for knockoff power and fdr estimation
q = 0.2; % target fdr

sigma = 1; % noise variance

% results
power_list_jlt = zeros(1,length(list_len));
fdr_list_jlt = zeros(1,length(list_len));
power_list_ag = zeros(3,length(list_len));
fdr_list_ag = zeros(3,length(list_len));


% Experiment
for i = 1:list_len
    [mean_power_jlt,fdr_jlt,mean_power_ag,fdr_ag] = DPknockoffs_compare(n_list(i),p_list,s0_list,mu_list,r_list,eps,del,b,n_iter,lambda,q);
    power_list_jlt(i) = mean_power_jlt;
    fdr_list_jlt(i) = fdr_jlt;
    power_list_ag(:,i) = mean_power_ag;
    fdr_list_ag(:,i) = fdr_ag;
end

save("DPknockoffs_expm_23.mat")
