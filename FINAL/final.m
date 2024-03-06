function [x, obj] = final(m, n, s, p)
    % Generates and solves a random non linear programming problem,
    % with the objective and constraints as specified in problem 1

    % Inputs:
    %   m,n - The dimensions of the problem
    %   s - The seed for the random number generator
    %   p - The maximum number of iterations allowed for the interior point algorithm
    
    % Outputs:
    %   x - The solution vector
    %   obj - The optimal objective function value

    % Nonconvexity parameter
    gamma = 0;
    
    % Generate a random linear program
    [c, Q0, D0, D, a] = RandomNonLinearProgram(m, n, s);
    
    % Solve the nonlinear program using the BCLMethod
    x = BCLMethod(c, Q0, D0, D, a, gamma, p);

    % Compute the value of the objective at x
    obj = c' * x + x' * Q0 * x + gamma * x' * D0 * x;
end
