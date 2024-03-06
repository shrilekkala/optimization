function x_cg = projected_CG(x_c, l, u, G, c)
    % Projected Conjugate Gradient Method (Algorithm 16.2)

    % INPUTS:
    % x_c - Current estimate of the solution
    % l - Lower bounds for x
    % u - Upper bounds for x
    % G - Approximate Hessian matrix (or true Hessian) at the current estimate
    % c - Gradient of the objective function at the current estimate

    % Initialization
    n = size(x_c, 1);   % Problem dimension
    I = eye(n);         % Identity matrix of the same size as G
    x_cg = x_c;         % Initialize x_cg with the current estimate
    r = G * x_cg + c;   % Residual vector
    
    % Determine active constraints and compute g
    active_constraint_idx = find(abs(x_cg - l) < 1e-3 | abs(x_cg - u) < 1e-3);
    A = I(active_constraint_idx, :);   

    % Solve for g via the augmented system approach (16.32)
    v = (A * A') \ (A * r);             
    g = r - A' * v;                     
    
    d = -g;
    cg_iter = 0;
    
    % CG iterations
    while r' * g > 1e-5 

        % Check for negative curvature
        if d' * G * d <= 0
            break;
        end
        
        % Step size calculation
        alpha = (r' * g) / (d' * G * d);
        x_prev = x_cg;
        x_cg = x_cg + alpha * d; % Update the estimate
        
        % Convergence check
        if norm(x_prev - x_cg) < 1e-10
            break;
        end
        
        % Enforce the bounds
        if any(x_cg < l) || any(x_cg > u)
            x_cg = x_prev; % Revert to the previous estimate if bounds violated
            break;
        end
        
        % Update residual and gradient
        r_plus = r + alpha * G * d;

        % Solve for g_plus via the augmented system approach (16.32)
        v = (A * A') \ (A * r_plus);
        g_plus = r_plus - A' * v;
        
        % Update search direction
        beta = (r_plus' * g_plus) / (r' * g);
        d = -g_plus + beta * d;
        g = g_plus;
        r = r_plus;
        
        % Update iteration counter
        cg_iter = cg_iter + 1;
    end
end
