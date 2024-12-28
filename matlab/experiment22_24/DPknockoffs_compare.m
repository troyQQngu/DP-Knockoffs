% DP Knockoffs Experiment kernel
function [mean_power_jlt,fdr_jlt,mean_power_ag,fdr_ag] = DPknockoffs_compare(n,p,s0,mu,r,eps,del,b,n_iter,lambda,q)

% nonparameters
B = sqrt(b^2*p*2); % bound on row norm of X
w_sq = B^2*(1+(2/eps+2/log(4/del))*(sqrt(2*r*log(4/del))+log(4/del)));
w = sqrt(w_sq);


% data initialization

S = sort(randperm(p,s0)); % set of true parameters
theta = zeros(p,1); % parameter vector
theta(S) = mu; % magnitude of the non zero parameters

% knockoffs iteration
% record jlt privatization results
fdps_jlt = zeros(1,n_iter);
powers_jlt = zeros(1,n_iter);

% record ag privatization results
% knockoffs iteration
fdps_ag = zeros(1,n_iter);
powers_ag = zeros(1,n_iter);
convs = zeros(1,n_iter);
mineigs = zeros(1,n_iter);


X = rand(n,p)*2*b-b; % design matrix


for i = 1:n_iter
    tic
    xi = randn(n,1);
    y = X*theta+xi; % response
    X_ko = rand(n,p)*2*b-b; % knockoffs

    std = B^2*sqrt(2*log(1.25/del))/eps;
    %std = 0;
    [XtX,Xty] = Analyzegauss(X,X_ko,y,std); % Gaussian Mechanism privatization
    [RX,RX_ko,Ry] = JLT(X,X_ko,y,r,w); % JLT privatization
    mineigs(i) = min(eig(XtX));
    

    %[W,theta_hat,lambda] = knockoffs.stats.lassoCoefDiffwlambda(RX,RX_ko,Ry,lambda_list); % compute statistics
    %grad = 1/n*[RX RX_ko]'*(Ry-[RX RX_ko]*theta_hat);
    %lambda = max(abs(grad)); % lambda calibration

    theta_hat_jlt = lasso_regression([RX RX_ko], Ry, lambda*n/r, 1e-8, 1000);
    [theta_hat_ag,convs(i)] = lasso_regression_adj(XtX,Xty,n,lambda,1e-8,1000);
    if convs(i) == 0
        theta_hat_ag = zeros(2*p,1);
    end

    W_jlt = abs(theta_hat_jlt(1:p))-abs(theta_hat_jlt((p+1):(2*p)));
    W_ag = abs(theta_hat_ag(1:p))-abs(theta_hat_ag((p+1):(2*p)));

    % for analysis puporses
    %********************************************************************************************************
    % theta_ap = [theta; zeros(p,1)];
    % theta_u = theta_hat_my+1/(gamma*n)*[RX RX_ko]'*(Ry-[RX RX_ko]*theta_hat_my)+w_sq/(n*gamma)*theta_hat_my;
    % tri = (M*[X X_ko]'*[X X_ko]-diag(ones(1,2*p)))*(theta_ap-theta_hat_my);
    % X_a = [X, X_ko; w*diag(ones(1,2*p));zeros(1,2*p)];
    % %R = [R1 R2 R3 R4];
    % %tri_til = 1/(n*gamma)*X_a'*(R'*R-diag(ones(1,n+2*p+1)))*[[X X_ko]*(theta_ap-theta_hat);-w*theta_hat;w];
    % terma = M*((R1*[X X_ko])'*R1*[X X_ko]-[X X_ko]'*[X X_ko])*(theta_ap-theta_hat_my);
    % termb = -M*w_sq*([R2 R3]'*[R2 R3]-eye(2*p))*theta_hat_my;
    % termc = M*w*[R2 R3]'*R1*[X X_ko]*(theta_ap-theta_hat_my);
    % termd = -M*w*[X X_ko]'*R1'*[R2 R3]*theta_hat_my;
    % terme = M*w*[X X_ko]'*R1'*R4;
    % termf = M*w_sq*[R2 R3]'*R4;
    % tri_til = terma+termb+termc+termd+terme+termf;
    % 
    % noise_term = 1/(n*gamma)*([X X_ko]'*R1'+w*[R2 R3]')*R1*xi;
    % rho = sqrt((n+r+1)/(n*r*gamma)+w_sq/(gamma^2*n*r));
    % fprintf('sup norm of tri = %.5f, sup norm of tri_til = %.5f, standard deviation of noise = %.5f, rho_n = %.5f\n', max(abs(tri)),max(abs(tri_til)),std(noise_term),rho)

    %********************************************************************************************************
    
    % data dependent threshold
    S_jlt = selection(W_jlt,q);
    S_ag = selection(W_ag,q);


    power_jlt = sum(theta(S_jlt) ~= 0)/s0; % power
    fdp_jlt = sum(theta(S_jlt) == 0)/max(1,length(S_jlt)); % FDP

    power_ag = sum(theta(S_ag) ~= 0)/s0; % power
    fdp_ag = sum(theta(S_ag) == 0)/max(1,length(S_ag)); % FDP

    powers_jlt(i) = power_jlt; % record power
    fdps_jlt(i) = fdp_jlt; % record fdp
    
    powers_ag(i) = power_ag; % record power
    fdps_ag(i) = fdp_ag; % record fdp

    time = toc; % time per iteration

    % supervision message
    fprintf('%s\n',repmat('-',1,50)) % horizontal line
    fprintf('%d / %d for n = %d, p = %d, s0 = %d, mu = %.4f, r = %d, epsilon = %.3f, w = %.1f, lambda = %.4f\n',i,n_iter,n,p,s0,mu,r,eps,w,lambda) % print parameters
    fprintf('True model: ') 
    fprintf('%d,',S)% true model
    fprintf('\nJLT DP knockoffs select: ') 
    fprintf('%d,',S_jlt) % knockoff selection
    fprintf('\nFDP = %.4f, Power = %.4f',fdp_jlt,power_jlt)
    fprintf('\nAG DP knockoffs select: ')
    fprintf('%d,',S_ag) % knockoff selection
    fprintf('\nFDP = %.4f, Power = %.4f\n',fdp_ag,power_ag)
    fprintf('Elapse time = %.4f\n',time)
end

mean_power_jlt = mean(powers_jlt);
fdr_jlt = mean(fdps_jlt);

mean_power_ag = zeros(3,1);
fdr_ag = zeros(3,1);

% selection only when PSD, average over everything
mean_power_ag(1) = sum(powers_ag(mineigs>0))/n_iter;
fdr_ag(1) = sum(fdps_ag(mineigs>0))/n_iter;
% average over everything
mean_power_ag(2) = mean(powers_ag);
fdr_ag(2) = mean(fdps_ag);
% average over convergent results only
if sum(convs==1)==0
    mean_power_ag(3) = 0;
    fdr_ag(3) = 0;
else
    mean_power_ag(3) = mean(powers_ag(convs==1));
    fdr_ag(3) = mean(fdps_ag(convs==1));
end

% supervision message
    fprintf('%s\n',repmat('*',1,50)) % horizontal stars
    fprintf('For n = %d, p = %d, s0 = %d, mu = %.4f, r = %d, epsilon = %.3f, w = %.1f, lambda = %.4f\n',n,p,s0,mu,r,eps,w,lambda) % print parameters
    fprintf('\nJLT: FDR = %.4f, Power = %.4f\n',fdr_jlt,mean_power_jlt)
    fprintf('AG: FDR = %.4f, Power = %.4f\n',[fdr_ag';mean_power_ag'])
    fprintf('%s\n',repmat('*',1,50)) % horizontal stars
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

function [XtX,Xty] = Analyzegauss(X,X_ko,y,std)
    [n,p] = size(X);
    A = [X,X_ko,y];

    noise = zeros(2*p+1);
    for i = 1:2*p
        noise(i,i+1:end) = std*randn(1,2*p+1-i);
    end
    noise = noise'+noise + diag(std*randn(1,2*p+1));

    AtA = A'*A+noise;
    XtX = AtA(1:2*p,1:2*p);
    Xty = AtA(1:2*p,2*p+1);
end



