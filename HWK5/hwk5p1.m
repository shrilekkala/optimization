function [x] = hwk5p1(x_in, tol, iter)
    % Finds a minimum of the n-dimensional Rosenbrock function 
    % using the damped BFGS method with backtracking line search

    % x_in: Initial point
    % iter: Maximum number of iterations
    % tol: Stopping tolerance
    
    % Use the Trust Region method for the Rosenbrock function
    [x_star, ~, ~] = damped_BFGS(@rosenbrocknfgH, x_in, tol, iter);
    
    x = x_star;
end
