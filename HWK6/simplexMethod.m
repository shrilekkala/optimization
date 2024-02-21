function [x, obj, isDegenerate, isUnbounded, n_iter] = simplexMethod(A, b, c, p)
    % Implements the simplex algorithm (13.1) for solving linear programming problems.
    
    % Inputs:
    %   A - Coefficient matrix of the constraints
    %   b - Right-hand side vector of the constraints
    %   c - Coefficients vector for the objective function
    %   p - Maximum number of iterations allowed

    % Outputs:
    %   x - Optimal variable values for the decision variables
    %   obj - Optimal objective function value
    %   isDegenerate - Boolean flag for degeneracy in the solution
    %   isUnbounded - Boolean flag for unboundedness in the solution
    %   n_iter - Number of iterations performed by the algorithm

    % Initialize variables
    isDegenerate = false;
    isUnbounded = false;
    n_iter = 0;
    m = size(A, 1); 
    n = size(A, 2);
    
    % Create initial basis B
    B_set = (n-m+1):n; % index set B
    N_set = 1:(n-m);   % index set N

    x = zeros(n, 1); % Initialize x as a vector of zeros with length n
    
    % Main simplex algorithm loop
    while n_iter < p
        n_iter = n_iter + 1;

        % construct matrices Basic and Nonbasic matrices
        B = A(:, B_set);
        N = A(:, N_set);
    
        % primal variable
        xB = B \ b;
        % Starting point
        if n_iter == 1
            xB = [ones(m, 1)];
        end

        % dual variable (lambda)
        lambda = B' \ c(B_set);
        
        % pricing (reduced costs for non-basic variables)
        sN = c(N_set) - N' * lambda;
        
        % check for optimality
        if all(sN >= 0)
            break; % optimal point found
        end
        
        % Determine entering variable (most negative reduced cost)
        [~, enteringIndex] = min(sN);
        q = N_set(enteringIndex);
        
        % Compute direction vector d
        d = B \ A(:, q);
        
        % Check for unboundedness
        if all(d <= 0)
            isUnbounded = true;
            disp("Problem is unbounded")
            break;
        end
        
        % Determine leaving variable
        ratios = xB ./ d;
        ratios(d <= 0) = inf;  % Ensure non-negativity
        [minRatio, leavingIndex] = min(ratios);

        % Check for exact ties in the minimum ratio test for degeneracy
        if sum(ratios == minRatio) > 1
            isDegenerate = true;
            disp("Degenerate case encountered")
            break;
        end

        x_q = minRatio;
        
        % Update basis and solution
        xB = xB - d * x_q;
        xB(leavingIndex) = x_q;
        
        % Update index sets
        N_set(enteringIndex) = B_set(leavingIndex);
        B_set(leavingIndex) = q;
    end
    
    if n_iter >= p
        disp("Max number of iterations reached before simplex method finished")
    end

    optimalBasis = B_set;
    
    % Compute the solution and optimal objectives from the optimal basis
    x(optimalBasis) = xB;
    obj = c' * x; 

end
