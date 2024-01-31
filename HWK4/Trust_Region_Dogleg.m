function [x_k, k, n_linear_systems, n_function_evals] = Trust_Region_Dogleg(fgH, x0, tol, max_iter)
    % fgH: Function that returns function value, gradient, and Hessian
    % x0: Starting point
    % tol: stopping criterion for gradient convergence
    % max_iter: maximum number of iterations for algorithm to run

    % Parameters for trust region algorithm
    delta_hat = 1;   % overall bound on step lengths
    eta = 1/8;       % trust region parameter

    % Initialize variables
    k = 0;                      % Iteration count
    delta_k = delta_hat / 2;    % Trust region radius
    x_k = x0;                   % Initial guess
    n_linear_systems = 0;
    n_function_evals = 0;

    % Algorithm 4.1, with dogleg approach
    while k < max_iter
        % Calculate f, g, and H at current point
        [f_k, g_k, B_k] = fgH(x_k);
        B_k = sparse(B_k);

        % Check for convergence in gradient norm
        if norm(g_k) < tol
            fprintf('Dogleg Converged in %d iterations.\n', k);
            return;
        end

        % Update function evaluation count
        n_function_evals = n_function_evals + 1;

        % Solve for pB using the dogleg method
        p_k = compute_dogleg_step(g_k, B_k, delta_k); 
        
        % Compute rho_k using the actual and predicted reductions (4.4)
        % rho_k close to 1 indicates a good model prediction
        rho_k = (f_k - fgH(x_k + p_k)) / (f_k - model_f(f_k, g_k, B_k, p_k));

        % Update linear system solve and function evaluation counts
        % (note we do not include evaluating the model function here)
        n_linear_systems = n_linear_systems + 1;
        n_function_evals = n_function_evals + 1; 
        
        % Update the trust region radius based on current performance
        if rho_k < 0.25
            % Radius is decreased to make the next steps more conservative
            delta_k = 0.25 * delta_k;
        else
            if rho_k > 0.75 && norm(p_k) == delta_k
                % Radius is increased to allow for potentially larger steps
                delta_k = min(2 * delta_k, delta_hat);
            else
                delta_k = delta_k;
            end
        end
        
        % Update the solution
        if rho_k > eta
            x_k = x_k + p_k;
        else
            x_k = x_k;
        end
        
        % Increment iteration count
        k = k + 1; 
    end
    warning('Max number of iterations reached without convergence.')
end


