function [beta, intercept] = lasso_cd(X, y, lambda, tol, max_iter)
    if nargin < 4 || isempty(tol),      tol      = 1e-6;    end
    if nargin < 5 || isempty(max_iter), max_iter = 1000;    end

    [n, p] = size(X);

    %— 1) Center y and X
    y_mean = mean(y);
    y_cent = y - y_mean;

    X_mean = mean(X,1);
    X_cent = X - X_mean;          % implicit broadcast in modern MATLAB

    %— 2) Compute column norms (for scaling), detect zero‐variance
    X_norm = sqrt(sum(X_cent.^2,1));  
    zero_var = (X_norm == 0);

    % avoid division by zero: give them norm=1 so Xs(:,j)=0
    X_norm_adj = X_norm;
    X_norm_adj(zero_var) = 1;

    % scale
    Xs = X_cent ./ X_norm_adj;

    %— 3) Initialize
    beta     = zeros(p,1);
    beta_old = beta;

    for iter = 1:max_iter
        for j = 1:p
            if zero_var(j)
                % no information in this column → keep beta(j)=0
                continue
            end

            % partial residual (leave j out)
            r_j = y_cent - Xs*beta + Xs(:,j)*beta(j);

            rho = Xs(:,j)' * r_j;

            % soft‐threshold update
            if rho < -lambda
                bj = (rho + lambda);
            elseif rho >  lambda
                bj = (rho - lambda);
            else
                bj = 0;
            end
            beta(j) = bj;  % denominator is 1 since ||Xs(:,j)||^2 = 1
        end

        % convergence check
        if norm(beta - beta_old, Inf) < tol
            break;
        end
        beta_old = beta;
    end

    %— 4) Rescale back to original units
    beta(zero_var) = 0;             % be explicit
    beta(~zero_var) = beta(~zero_var) ./ X_norm(~zero_var)';

    %— 5) Intercept
    intercept = y_mean - X_mean * beta;
end
