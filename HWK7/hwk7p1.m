function [x, obj] = hwk7p1(m, n, s, p)
    % Generates and solves a random linear programming problem.

    % Inputs:
    %   m,n - The dimensions of the problem
    %   s - The seed for the random number generator
    %   p - The maximum number of iterations allowed for the interior point algorithm
    
    % Outputs:
    %   x - The solution vector
    %   obj - The optimal objective function value
    
    % Generate a random linear program
    [A, b, c] = RandomLinearProgram(m, n, s);
    
    % Solve the linear program using the simplex method
    [x, obj] = longStepInteriorPoint(A, b, c, p);
end
