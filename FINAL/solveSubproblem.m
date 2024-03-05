% Solution to subproblem (17.50) from the bound-constrained Lagrangian method.
function [x_k_plus, success] = solveSubproblem(x_k, lambda, mu, c, Q0, D0, D, a, gamma, omega)
    % Set algorithm parameters
    alpha_init = 1;     % Initial step size for line search
    rho = 0.9;          % Step size reduction factor
    c1 = 1e-4;          % Constant for Armijo rule in line search
    maxIter = 100;      % Maximum number of iterations
    
    % Initialize Boolean success flag
    success = true;
    
    % Define augmented Lagrangian and its gradient and Hessian
    L_A = @(x) augmentedLagrangian(x, lambda, mu, c, Q0, D0, D, a, gamma);
    grad_L_A = @(x) gradLagrangian(x, lambda, mu, c, Q0, D0, D, a, gamma);
    H_L_A = @(x) computeHessian(Q0, D0, D, a, gamma, lambda, mu, x);

    % Iterative process to update the current point
    for iter = 1:maxIter

        % Find an approximate minimizer of the quadratic subproblem
        % % x_hat = solveQuadraticSubproblem(x_k, grad_L_A(x_k), H_L_A(x_k));
        x_hat = gradient_projection_method(x_k, grad_L_A(x_k), H_L_A(x_k));

        % Check for NaNs and exit if found
        if any(isnan(x_hat))
            fprintf("Warning: NaN encountered\n");
            success = false;
            x_k_plus = x_k;
            return;
        end

        % Calculate the search direction
        p_k = x_hat - x_k;

        % Check for convergence based on the projected gradient norm
        if norm(x_k - projectOntoBounds(grad_L_A(x_k), -1, 1), 'inf') < omega
            break; % Convergence achieved
        end

        % Perform line search with Armijo rule
        alpha = alpha_init;
        while true
            x_k_new = x_k + alpha * p_k;
            if L_A(x_k_new) <= L_A(x_k) + c1 * alpha * grad_L_A(x_k)' * p_k
                break; % Sufficient decrease achieved
            else
                alpha = rho * alpha; % Reduce step size
            end
            if alpha < 1e-8
                break; % Avoiding infinite loop due to too small alpha
            end
        end

        % Update the current estimate
        x_k = x_k_new;
    end

    % Output the final estimate after iteration
    x_k_plus = x_k;
end
