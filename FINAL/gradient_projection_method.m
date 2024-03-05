% Algorithm 16.5 (Gradient Projection Method for QP)
function x_star = gradient_projection_method(x_k, g, H)
    % INPUTS:
    % x_k - The current point in the search space.
    % g - The gradient of the objective function at the initial point.
    % H - Positive definite matrix, an approximation to the Hessian.

    % Initializations
    n = size(x_k, 1);
    maxIter = 1000;
    
    % Set bounds for the variables
    l = -ones(n, 1);  % Lower bounds
    u = ones(n, 1);   % Upper bounds
    
    for k = 1:maxIter

        % Compute the gradient of the quadratic at the current point
        g_k = H * x_k + g;
        
        % Check if current point satisfies KKT conditions
        if checkKKT(x_k, g_k, l, u)
            x_star = x_k;  % Optimal point found
            return;
        end
        
        % Compute the Cauchy point (steepest descent direction)
        x_c = computeCauchyPoint(x_k, g_k, H, l, u);
        
        % Find an approximate solution x_plus using projected CG
        x_plus = projected_CG(x_c, l, u, H, g_k);
        
        % Compare objective function values at x_plus and x_c
        if q_k(x_plus, x_k, g_k, H) >= q_k(x_c, x_k, g_k, H)
            x_star = x_c;  % No improvement found
            return;
        end
        
        % Update the current point for the next iteration
        x_k = x_plus;
    end
    
    % Maximum number of iterations reached without satisfying KKT conditions
    x_star = x_k;
end

function isKKT = checkKKT(x, g, l, u)
    % Check KKT conditions for the bound-constrained problem
    % x - Current point
    % g - Gradient of the objective function at x
    % l - Lower bounds
    % u - Upper bounds

    tol = 1e-4; % Tolerance for the stopping criterion
    
    % Initialize the check as true
    isKKT = true;

    % Check conditions for each component
    for i = 1:length(x)
        if (x(i) > l(i) && x(i) < u(i) && abs(g(i)) > tol)
            isKKT = false; % Gradient should be zero in the interior
            break;
        elseif (x(i) == l(i) && g(i) < -tol)
            isKKT = false; % Gradient should be non-negative at the lower bound
            break;
        elseif (x(i) == u(i) && g(i) > tol)
            isKKT = false; % Gradient should be non-positive at the upper bound
            break;
        end
    end
end

function val = q_k(x, x_k, g_k, B_k)
    % Evaluate the quadratic model q_k at the point x
    val = (x-x_k)' * g_k + 0.5 * (x-x_k)' * B_k * (x-x_k);
end