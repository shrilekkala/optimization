function [x, obj] = longStepInteriorPoint(A, b, c, p)
    % Input parameters:
    % A, b, c: Coefficients of the linear program
    % p: Maximum number of iterations allowed

    % Outputs:
    % x - Optimal variable values for the decision variables
    % obj - Optimal objective function value

    % Initialize parameters
    [m, n] = size(A);
    e = ones(n, 1);
    x = ones(n, 1);           % initial x is the vector of ones
    lambda = zeros(m, 1);     % initial lambda is the zero vector
    s = c;                    % initial s is the vector c

    % Set convergence tolerance
    tol = 1e-6;

    % Set sigma_min and sigma_max
    sigma_min = 0.1;
    sigma_max = 0.9;

    % Calculate mu and gamma
    mu = (x' * s) / n;
    gamma = min(s ./ mu);
    
    for k = 1:p
        
        % Choose sigma_k
        sigma_k = sigma_max;

        % Calculate residuals
        r_c = A' * lambda + s - c;
        r_b = A * x - b;

        % Compute the duality measure
        mu = (x' * s) / n;
        
        % Construct the KKT matrix
        S = diag(s);           % diagonal matrix of s
        X = diag(x);           % diagonal matrix of x
        KKT = [zeros(n, n), A'         ,    eye(n)      ;
               A          , zeros(m, m),    zeros(m, n) ;
               S          , zeros(n, m),    X          ];
        
        % Construct the right-hand side vector
        rhs = [-r_c;
               -r_b;
               -X * S * e + sigma_k * mu * e];

        % (14.10) Solve for search directions from KKT * Δ = rhs
        deltas = KKT \ rhs;
        delta_x = deltas(1:n);
        delta_lambda = deltas(n+1:n+m);
        delta_s = deltas(n+m+1:end);
        
        % Compute the step length α
        alpha = lineSearchAlgorithm(x, s, delta_x, delta_s, gamma, mu);
        
        % Update x, lambda, and s
        x = x + alpha * delta_x;
        lambda = lambda + alpha * delta_lambda;
        s = s + alpha * delta_s;

        % Stopping criterion for convergence
        if norm(r_c, Inf) < tol && norm(r_b, Inf) < tol && norm(mu) < tol;
            fprintf('Interior point converged in %d iterations.\n', k);
            obj = c' * x;
            return;
        end
    end

    % Compute the objective
    obj = c' * x

    % If loop has exited, then maximum number of iterations have been reached
    warning('Max number of iterations reached without convergence.');
    return;
end