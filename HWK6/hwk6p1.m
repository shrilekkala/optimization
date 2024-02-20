function [x, obj] = hwk6p1(m, n, s, p)
    % Generates and solves a random linear programming problem.

    % Inputs:
    %   m,n - The dimensions of the problem
    %   s - The seed for the random number generator
    %   p - The maximum number of iterations for the simplex method
    
    % Outputs:
    %   x - The solution vector
    %   obj - The optimal objective function value
    
    % Generate a random linear program
    [A, b, c] = RandomLinearProgram(m, n, s);
    
    % Solve the linear program using the simplex method
    [x, obj, ~, ~, ~] = simplexMethod(A, b, c, p);
end
