function [x] = hwk3p3(x0, niter, eps)
    % Finds a minimum of the n-dimensional Rosenbrock function using NMHM

    % x0: Initial point
    % niter: Maximum number of iterations
    % eps: Gradient precision for stopping criterion
    
    % Use the NMHM method for the Rosenbrock function
    % Note in the NMHM function there is a double test for exiting the loop
    x = NMHM(@rosenbrocknfgH, x0, eps, niter);
end
