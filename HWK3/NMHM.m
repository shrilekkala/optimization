function x = NMHM(fgH, x0, tol, max_iter)
    % fgH: Function that returns function value, gradient, and Hessian
    % Starting point
    x = x0;
    
    % Parameters for backtracking line search
    alpha_bar = 1;  % Initial step size
    rho = 0.5;      % Step size reduction factor
    c = 1e-4;       % Armijo rule constant
    
    % Parameters for Hessian modification
    beta = 1e-3;    % Hessian modification factor
    
    % Algorithm 3.2 (Line Search Newton with Modification)
    for k = 1:max_iter

        % Calculate gradient and Hessian at current point
        [f_val, g, H] = fgH(x);
        
        % Modify the Hessian matrix using Algorithm 3.3
        Bk = modify_hessian_with_cholesky(H, beta);
        
        % Solve Bk * pk = - Grad(f(xk)) for pk
        pk = - Bk \ g;
        
        % Find step size using backtracking Algorithm 3.1
        alpha_k = backtracking_line_search(fgH, x, pk, alpha_bar, rho, c);
        
        % Update the point
        x = x + alpha_k * pk;
        
        % Check for convergence
        if norm(g) < tol
            fprintf('NMHM Converged in %d iterations.\n', k);
            return;
        end
    end
    warning('Max number of iterations reached without convergence.')
end

function Bk = modify_hessian_with_cholesky(H, beta)
    % Algorithm 3.3
    if min(diag(H)) > 0
        tau = 0;
    else
        tau = -min(diag(H)) + beta;
    end
    
    k = 0;
    while true
        % apply Cholesky algorithm
        [L, p] = chol(H + tau * eye(size(H)), 'lower');
        % p = 0 signifies that the factorization was successful
        if p == 0
            Bk = L * L';
            return;
        else
            tau = max(2 * tau, beta);
        end
        k = k + 1;
        % exit condition
        if k > 100
            error('Hessian could not be modified to be positive definite.');
        end
    end
end

function alpha_k = backtracking_line_search(fgH, x, pk, alpha_bar, rho, c)
    % Algorithm 3.1
    alpha = alpha_bar;

    [f_next, ~, ~] = fgH(x + alpha * pk);
    [f_current, grad_current, ~] = fgH(x);

    while f_next > f_current + c * alpha * grad_current' * pk
        alpha = rho * alpha; % Contract the step size
        % update f_next using the new alpha:
        [f_next, ~, ~] = fgH(x + alpha * pk);
    end
    alpha_k = alpha;
end
