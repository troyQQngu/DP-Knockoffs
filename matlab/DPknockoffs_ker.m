% DP Knockoffs Experiment kernel
function [mean_power,fdr,power_th,fdr_th,fdr_est_th] = DPknockoffs_ker(n,p,s0,mu,r,eps,del,b,n_iter,T,lambda,gamma,c0)
% nonparameters
M = 1/(gamma*n);
B = sqrt(b^2*p); % bound on row norm of X
w_sq = B^2*(1+(2/eps+2/log(4/del))*(sqrt(2*r*log(4/del))+log(4/del)));
%w_sq = p*sqrt(r)/eps;
w = sqrt(w_sq);


% data initialization

S = sort(randperm(p,s0)); % set of true parameters;
theta = zeros(p,1); % parameter vector
theta(S) = mu; % magnitude of the non zero parameters


%% knockoffs iteration
fdps = zeros(1,n_iter);
powers = zeros(1,n_iter);
lambdas = zeros(1,n_iter);

X = rand(n,p)*2*b-b; % design matrix

for i = 1:n_iter
    tic
    xi = randn(n,1);
    y = X*theta+xi; % response
    X_ko = rand(n,p)*2*b-b; % knockoffs
    [RX,RX_ko,Ry,] = JLT(X,X_ko,y,r,w); % JLT privatization

    %[W,theta_hat,lambda] = knockoffs.stats.lassoCoefDiffwlambda(RX,RX_ko,Ry,lambda_list); % compute statistics
    %grad = 1/n*[RX RX_ko]'*(Ry-[RX RX_ko]*theta_hat);
    %lambda = max(abs(grad)); % lambda calibration

    theta_hat_my = lasso_regression([RX RX_ko], Ry, lambda*n/r, 1e-10, 1000); 
    W = abs(theta_hat_my(1:p))-abs(theta_hat_my((p+1):(2*p)));

    % for analysis puporses
    %********************************************************************************************************
    % theta_ap = [theta; zeros(p,1)];
    % theta_u = theta_hat+1/(gamma*n)*[RX RX_ko]'*(Ry-[RX RX_ko]*theta_hat)+w_sq/(n*gamma)*theta_hat;
    % tri = (M*[X X_ko]'*[X X_ko]-diag(ones(1,2*p)))*(theta_ap-theta_hat);
    % X_a = [X, X_ko; w*diag(ones(1,2*p));zeros(1,2*p)];
    % %R = [R1 R2 R3 R4];
    % %tri_til = 1/(n*gamma)*X_a'*(R'*R-diag(ones(1,n+2*p+1)))*[[X X_ko]*(theta_ap-theta_hat);-w*theta_hat;w];
    % terma = M*((R1*[X X_ko])'*R1*[X X_ko]-[X X_ko]'*[X X_ko])*(theta_ap-theta_hat);
    % termb = -M*w_sq*([R2 R3]'*[R2 R3]-eye(2*p))*theta_hat;
    % termc = M*w*[R2 R3]'*R1*[X X_ko]*(theta_ap-theta_hat);
    % termd = -M*w*[X X_ko]'*R1'*[R2 R3]*theta_hat;
    % terme = M*w*[X X_ko]'*R1'*R4;
    % termf = M*w_sq*[R2 R3]'*R4;
    % tri_til = terma+termb+termc+termd+terme+termf;
    % 
    % noise_term = 1/(n*gamma)*([X X_ko]'*R1'+w*[R2 R3]')*R1*xi;
    % rho = sqrt((n+r+1)/(n*r*gamma)+w_sq/(gamma^2*n*r));
    % fprintf('sup norm of tri = %.5f, sup norm of tri_til = %.5f, standard deviation of noise = %.5f, rho_n = %.5f\n', max(abs(tri)),max(abs(tri_til)),std(noise_term),rho)

    %********************************************************************************************************
    
    S_hat = find(W>=T); % selection
    power = sum(theta(S_hat) ~= 0)/s0; % power
    fdp = sum(theta(S_hat) == 0)/max(1,length(S_hat)); % FDP

    powers(i) = power; % record power
    fdps(i) = fdp; % record fdp
    lambdas(i) = lambda; % record lambda
    time = toc; % time per iteration

    % supervision message
    fprintf('%s\n',repmat('-',1,50)) % horizontal line
    fprintf('%d / %d for n = %d, p = %d, s0 = %d, mu = %.4f, r = %d, epsilon = %.1f, w = %.1f, T = %.4f\n',i,n_iter,n,p,s0,mu,r,eps,w,T) % print parameters
    fprintf('DP knockoffs select: ') 
    fprintf('%d,',S_hat) % knockoff selection
    fprintf('\nTrue model: ') 
    fprintf('%d,',S)% true model
    fprintf('\nFDP = %.4f, Power = %.4f, lambda = %.4f\n',fdp,power,lambda)
    fprintf('Elapse time = %.4f\n',time)
end

mean_power = mean(powers);
fdr = mean(fdps);
mean_lambda = mean(lambdas);

% theoretical part
eta = @(x) 1/(1+w_sq/(gamma*n))*sign(x).*max((abs(x)-mean_lambda/gamma),0);
rho = sqrt((n+r+1)/(n*r*gamma)+w_sq/(gamma^2*n*r));

% monte carlo 
P1 = 0;
P2 = 0;
P3 = 0;
mc_iter = 1000;
for i = 1: mc_iter
    Z = randn(1,6);
    if abs(eta(mu+rho*Z(1)))-abs(eta(rho*Z(2)))<=-T
        P1 = P1+1;
    end
    if abs(eta(mu+rho*Z(3)))-abs(eta(rho*Z(4)))>=T
        P2 = P2+1;
    end
    if abs(eta(rho*Z(5)))-abs(eta(rho*Z(6)))<=-T
        P3 = P3+1;
    end
end
P1 = P1/mc_iter;
P2 = P2/mc_iter;
P3 = P3/mc_iter;
fdr_est_th = (c0*P1+(1-c0)*P3)/(c0*P2+(1-c0)*P3);
fdr_th = (1-c0)*P3/(c0*P2+(1-c0)*P3);
power_th = P2;

fprintf('%s\n',repmat('*',1,50)) % horizontal line 
fprintf('For n = %d, p = %d, s0 = %d, mu = %.4f, r = %d, epsilon = %.1f, w = %.1f, T = %.4f\n',n,p,s0,mu,r,eps,w,T)
fprintf('FDR = %.2f, Average power = %.2f, Average lambda_min = %.4f\n', fdr,mean_power,mean_lambda)
fprintf('Theoretical FDR = %.2f, Theoretical power = %.2f\n', fdr_th,power_th)
fprintf('%s\n',repmat('*',1,50)) % horizontal line

end

function [RX,RX_ko,Ry] = JLT(X,X_ko,y,r,w)
    [n,p] = size(X);

    % privatization
    RX = zeros(r,p);
    RX_ko = zeros(r,p);
    Ry = zeros(r,1);
    patchsize = 500;
    n_iter = floor(r/patchsize);
    
    for i = 1:n_iter
        R1 = randn(patchsize,n)/sqrt(r);
        R2 = randn(patchsize,p)/sqrt(r);
        R3 = randn(patchsize,p)/sqrt(r);
        R4 = randn(patchsize,1)/sqrt(r);
        RX((i-1)*patchsize+1:i*patchsize,:) = R1*X+w*R2;
        RX_ko((i-1)*patchsize+1:i*patchsize,:) = R1*X_ko+w*R3;
        Ry((i-1)*patchsize+1:i*patchsize) = R1*y+w*R4;
    end

    if i*patchsize < r
        lastpatch = r-i*patchsize;
        R1 = randn(lastpatch,n)/sqrt(r);
        R2 = randn(lastpatch,p)/sqrt(r);
        R3 = randn(lastpatch,p)/sqrt(r);
        R4 = randn(lastpatch,1)/sqrt(r);
        RX(i*patchsize+1:r,:) = R1*X+w*R2;
        RX_ko(i*patchsize+1:r,:) = R1*X_ko+w*R3;
        Ry(i*patchsize+1:r) = R1*y+w*R4;
    end
end




