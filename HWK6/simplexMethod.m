function [x, obj, isDegenerate, isUnbounded, n_iter] = simplexMethod(A, b, c, p)
    
    % Initialize variables
    isDegenerate = false;
    isUnbounded = false;
    n_iter = 0;
    m = size(A, 1); 
    n = size(A, 2);
    
    % Create initial basis B
    % Assuming that A is in standard form and the last m columns can be an identity matrix
    B_set = (n-m+1):n; % index set B
    N_set = 1:(n-m);   % index set N
    
    % Main simplex algorithm loop
    while n_iter < p
        n_iter = n_iter + 1;

        % construct matrices Basic and Nonbasic matrices
        B = A(:, B_set);
        N = A(:, N_set);
    
        % primal variable
        xB = B \ b;
        
        % Check for feasibility
        if any(xB < 0)
            error('Initial basis is not feasible');
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
        x_q = minRatio;
        
        % Check for degeneracy
        if x_q == inf
            isDegenerate = true;
            disp("Degenerate case encountered")
            break;
        end
        
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
    x = xB(optimalBasis);
    obj = c(1:m)' * x;
end
