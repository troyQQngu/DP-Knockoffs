% simulate
n=150; p=20;
X = randn(n,p);
btrue = [2; -3; zeros(p-2,1)];
y = X*btrue + randn(n,1)*0.3;

% 5-fold CV with MATLAB’s default λ-grid
[lamstar, bhat, b0hat, mse_mean, lamgrid] = lasso_cv(X, y);

fprintf('Chosen λ = %.4g\n', lamstar);
disp('Coefficients:'); disp(bhat);
