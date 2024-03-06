function isSOSC = checkSOSC(x, lambda, Q0, D0, D, a, gamma)
    % Compute the Hessian of the Lagrangian
    H_L = lagrangianHessian(x, lambda, Q0, D0, D, a, gamma);
    
    % Identify active constraints
    [activeEquality, activeBounds] = identifyActiveConstraints(x, a, D, gamma);
    
    % Project the Hessian onto the tangent space of the active constraints
    A = [activeEquality; activeBounds];

    % Otherwise, project Hessian onto the tangent space
    P = eye(length(x)) - A' * ((A * A') \ A);
    reducedHessian = P' * H_L * P;
    
    % Check if the reduced Hessian is positive definite
    eigenvalues = eig(reducedHessian);
    isSOSC = all(eigenvalues > 0);
end

function [activeEquality, activeBounds] = identifyActiveConstraints(x, a, D, gamma)
    m = size(a, 1);  % Number of equality constraints
    n = length(x);   % Number of variables
    
    % Tolerance for identifying active bound constraints
    tol = 1e-6;
    
    % Equality constraints are always active, so we calculate their gradients
    activeEquality = zeros(m, n);
    for i = 1:m
        activeEquality(i, :) = a(i, :) + 2 * gamma * (D{i} * x)';
    end

    % Identify active bound constraints
    activeBounds = [];
    for i = 1:n
        % Bounds are -1 < x < 1
        if abs(x(i) + 1) <= tol || abs(x(i) - 1) <= tol

            % If x(i) is within tol of the bound, it's an active constraint.
            % Create a unit vector with 1 in the position of the active bound
            unitVec = zeros(n, 1);
            unitVec(i) = 1;
            activeBounds = [activeBounds; unitVec'];
        end
    end
end


function H_L = lagrangianHessian(x, lambda, Q0, D0, D, a, gamma)
    % Dimensions
    m = size(a, 1);
    n = length(x);

    % Hessian of the objective function
    H_f = 2 * Q0 + D0;

    % Hessian for the augmented terms related to equality constraints
    H_c = zeros(n, n);
    for i = 1:m
        % Hessian of the equality constraint
        H_c = H_c + lambda(i) * (2 * gamma * D{i});
    end

    % Combine them to form the Hessian of the Lagrangian
    H_L = H_f + H_c;
end


