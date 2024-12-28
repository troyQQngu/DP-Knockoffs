function [beta,conv] = lasso_regression_adj(XtX, XtY,n, lambda, tol, max_iter)
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
    p = size(XtX,1);
    beta = zeros(p, 1);  % Start with zero coefficients
    beta_old = beta;     % Copy to check for convergence
    iter = 0;
    conv=1;
   
    % Precompute some values
    X_norm = diag(XtX);  % Column-wise squared sum of X

    % Main coordinate descent loop
    while iter < max_iter
        iter = iter + 1;
       
        % Coordinate descent
        for j = 1:p
            % Update beta_j using soft-thresholding
            rho = 1/n*(XtY(j) - XtX(j,:) * beta + XtX(j,j) * beta(j));  % Correlation term

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
        conv = 0;
    end
end