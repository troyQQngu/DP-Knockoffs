function [bestLambda, beta_best, intercept_best, mse_mean, lambda_seq] = lasso_cv(...
        X, y, lambda_seq, K, tol, max_iter)
% lasso_cv   K‐fold CV for our LASSO (coordinate descent), with
%            built-in default λ‐grid if you pass lambda_seq = [] or omit it.
%
% [bestLambda, β, b0, mse_mean, lambda_seq] = ...
%     lasso_cv(X, y)
%
% [bestLambda, β, b0, mse_mean, lambda_seq] = ...
%     lasso_cv(X, y, lambda_seq, K, tol, max_iter)
%
    if nargin < 4 || isempty(K),        K        = 5;      end
    if nargin < 5 || isempty(tol),      tol      = 1e-6;   end
    if nargin < 6 || isempty(max_iter), max_iter = 1000;   end

    [n, ~] = size(X);

    %— generate default λ‐grid if needed
    if nargin < 3 || isempty(lambda_seq)
        lambda_seq = defaultLambdaGrid(X, y, 100, 1e-4);
    end

    L        = numel(lambda_seq);
    mse      = zeros(K, L);
    idx      = randperm(n);
    foldSize = floor(n/K);

    %— K‐fold CV
    for i = 1:K
        if i < K
            valIdx = idx((i-1)*foldSize + (1:foldSize));
        else
            valIdx = idx((i-1)*foldSize + 1 : n);
        end
        trainMask = true(n,1);
        trainMask(valIdx) = false;

        Xtr = X(trainMask,:);  ytr = y(trainMask);
        Xva = X(~trainMask,:); yva = y(~trainMask);

        for j = 1:L
            lam = lambda_seq(j);
            [b_cv, b0_cv] = lasso_cd(Xtr, ytr, lam, tol, max_iter);
            ypred = Xva*b_cv + b0_cv;
            mse(i,j) = mean((yva - ypred).^2);
            fprintf('%.1f/%.1f\n',i,j)
        end
    end

    mse_mean      = mean(mse,1);
    [~, idxMin]   = min(mse_mean);
    bestLambda    = lambda_seq(idxMin);

    %— final fit on full data
    [beta_best, intercept_best] = lasso_cd(X, y, bestLambda, tol, max_iter);
end
