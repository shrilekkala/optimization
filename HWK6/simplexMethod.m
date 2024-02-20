function [optimalBasis, isDegenerate, optimalValue, isUnbounded, n_iter] = simplexMethod(A, b, c, p)
    
    % Initialize variables
    isDegenerate = false;
    isUnbounded = false;
    n_iter = 0;
    m = size(A, 1); n = size(A, 2);
    
    % Create initial basis B.
    B_set = m+1:n; % index set B
    N_set = 1:m;   % index set N
    
    
    % Main simplex algorithm loop
    while n_iter < p
        B = A(:, B_set); % matrix B
        N = A(:, N_set); % matrix N
    
        xB = B \ b;  % Initial basic feasible solution
        xB

        B_set
        N_set
        n_iter = n_iter + 1;
        
        lambda = B' \ c(B_set);
        lambda
        
        % Pricing
        sN = c(N_set) - N' * lambda;
        c(N_set)
        A
        N
        N' * lambda
        sN
        
        % Check for optimality
        if all(sN >= 0)
            optimalValue = c(B_set)' * xB;
            optimalBasis = B_set;
            return;
        end
        
        % Determine entering variable (most negative reduced cost)
        [minReducedCost, enteringIndex] = min(sN);
        if minReducedCost >= 0  % If no negative reduced cost, solution is optimal
            break;
        end
        
        q = N_set(enteringIndex);
        q
        
        % Compute direction vector d
        d = B \ A(:, q);
        d
        
        % Check for unboundedness
        if all(d <= 0)
            isUnbounded = true;
            break;
        end
        
        % Determine leaving variable
        ratios = xB ./ d;
        ratios(d <= 0) = inf;  % Ensure non-negativity of d_i
        [minRatio, leavingIndex] = min(ratios);
        
        x_q = minRatio
        p = leavingIndex
        

        % Check for degeneracy
        if x_q == 1
            isDegenerate = true;
            break;
        end
        
        % Update basis and solution
        xB = xB - d * x_q;
        xN(p) = x_q;

        % Update index sets
        N_set(enteringIndex) = B_set(p);
        B_set(p) = q;
        
    end
    
    % If the loop exits without returning, the maximum number of iterations was reached
    optimalValue = c(B_set)' * xB;
    optimalBasis = B_set;
end
