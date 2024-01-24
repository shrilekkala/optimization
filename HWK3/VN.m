function x = VN(fgH, x0, tol, max_iter)
    % fgh: Function that returns function value, gradient, and Hessian
    % Initialize starting point
    x = x0;

    for k = 1:max_iter
        % Calculate gradient and Hessian at current point
        [~, g, H] = fgH(x);
        
        % Check if Hessian is singular
        if det(H) == 0
            error('Hessian is singular');
        end

        % Update the solution
        x = x - H \ g;
        
        % Check for convergence
        if norm(g) < tol
            fprintf('VN Converged in %d iterations.\n', k);
            return;
        end
    end
    
    warning('Max number of iterations reached without convergence.');
end
