function beta = lasso_regression(X, Y, lambda, tol, max_iter)
    % lasso_regression: Perform LASSO regression using coordinate descent
    %
    % Inputs:
    %   X - n x p matrix of predictors
    %   Y - n x 1 vector of responses
    %   lambda - regularization parameter
    %   tol - tolerance for convergence
    %   max_iter - maximum number of iterations
    %
    % Output:
    %   beta - p x 1 vector of LASSO coefficients

    % Initialize parameters
    [n, p] = size(X);
    beta = zeros(p, 1);  % Start with zero coefficients
    beta_old = beta;     % Copy to check for convergence
    iter = 0;
   
    % Precompute some values
    X_norm = sum(X.^2);  % Column-wise squared sum of X

    % Main coordinate descent loop
    while iter < max_iter
        iter = iter + 1;
       
        % Coordinate descent
        for j = 1:p
            % Calculate the partial residual, excluding the current beta_j
            r_j = Y - X * beta + X(:, j) * beta(j);
           
            % Update beta_j using soft-thresholding
            rho = 1/n*X(:, j)' * r_j;  % Correlation term
            if rho < -lambda
                beta(j) = (rho + lambda) / (X_norm(j)/n);
            elseif rho > lambda
                beta(j) = (rho - lambda) / (X_norm(j)/n);
            else
                beta(j) = 0;
            end
        end
       
        % Check for convergence
        if norm(beta - beta_old, 2) < tol
            disp(['Converged in ' num2str(iter) ' iterations']);
            break;
        end
        beta_old = beta;
    end
   
    if iter == max_iter
        disp('Reached maximum iterations without full convergence');
    end
end