function [x] = hwk4p1(x0, niter, eps)
    % Finds a minimum of the n-dimensional Rosenbrock function 
    % using the Trust Region Dogleg method

    % x0: Initial point
    % niter: Maximum number of iterations
    % eps: Gradient precision for stopping criterion
    
    % Use the Trust Region method for the Rosenbrock function
    [x_star, n_iters, n_systems, n_evals] = Trust_Region_Dogleg(@rosenbrocknfgH, x0, eps, niter);

    x = x_star;
end
