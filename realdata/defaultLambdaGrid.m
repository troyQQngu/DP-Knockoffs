function lambda_seq = defaultLambdaGrid(X, y, numLambda, lambdaRatio)
% defaultLambdaGrid  Build a default λ‐grid like MATLAB’s lasso does
%
%   lambda_seq = defaultLambdaGrid(X, y)
%   lambda_seq = defaultLambdaGrid(X, y, numLambda)
%   lambda_seq = defaultLambdaGrid(X, y, numLambda, lambdaRatio)
%
% Inputs:
%   X            n×p data matrix
%   y            n×1 response
%   numLambda    number of λ’s (default 100)
%   lambdaRatio  min λ as fraction of max λ (default 1e-4)
%
% Output:
%   lambda_seq   numLambda×1 vector, geometric from λ_max→λ_max*lambdaRatio

    if nargin < 3 || isempty(numLambda),   numLambda   = 100;   end
    if nargin < 4 || isempty(lambdaRatio), lambdaRatio = 1e-4;  end

    % 1) center y
    y_cent = y - mean(y);

    % 2) center & scale X
    X_cent  = X - mean(X,1);
    X_norms = sqrt(sum(X_cent.^2,1));
    Xs      = X_cent ./ X_norms;        % each column has unit ℓ2‐norm

    % 3) compute λ_max = max_j |x_j' y|
    lambda_max = max(abs(Xs' * y_cent));

    % 4) λ_min
    lambda_min = lambda_max * lambdaRatio;

    % 5) geometric sequence
    lambda_seq = exp( linspace(log(lambda_max), log(lambda_min), numLambda) )';
end
