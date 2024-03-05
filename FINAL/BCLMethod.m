% Implements the Bound-Constrained Lagrangian (BCL) method for optimization
function [x, lambda, mu] = BCLMethod(c, Q0, D0, D, a, gamma, maxIter)

    % INPUTS:
    % c, Q0, D0 - vector and matrices of the objective function
    % a, D      - vector and matrices of the constraints
    % gamma     - nonconvecity paramter
    % maxIter - Maximum number of iterations

    % OUTPUTS:
    % x      - Solution vector
    % lambda - Lagrange multipliers for equality constraints
    % mu     - Penalty parameter

    % Set the bounds for the problem
    l = -1;
    u = 1;

    % Initialize convergence tolerances
    eta_star = 1e-4;
    omega_star = 1e-4;

    % Initialize variables
    n = length(c);      % Number of variables
    m = size(a, 1);     % Number of constraints
    x = zeros(n, 1);    % Start at the origin
    lambda = zeros(m, 1); % Start with zero multipliers

    % Initialize penalty parameter
    mu = 10;
    omega = 1 / mu;
    eta = 1 / (mu^0.1);

    % Iterative optimization process
    for k = 1:maxIter
        fprintf('Iteration: %d\n', k);
        obj_value = computeObjective(x, c, Q0, D0, gamma);
        fprintf('Objective value: %.4f\n', obj_value);

        % Solve the subproblem for the current iteration
        [x_new, success] = solveSubproblem(x, lambda, mu, c, Q0, D0, D, a, gamma, omega);

        % Exit if subproblem solution failed
        if ~success
            fprintf("Warning: Method terminated due to failure in solving the subproblem\n");
            return;
        end

        x = x_new; % Update solution

        % Compute constraint violations
        c_val = calculateC(x, a, D, gamma);
        proj_grad_norm = norm(x - projectOntoBounds(gradLagrangian(x, lambda, mu, c, Q0, D0, D, a, gamma), l, u));
    
        % Adjust penalty parameters based on constraint violations
        if norm(c_val) <= max(eta, eta_star)

            % Check for convergence
            if norm(c_val) <= eta_star && proj_grad_norm <= omega_star
                fprintf("Convergence achieved after %d iterations.\n", k);
                break; % Convergence achieved
            end
            
            % Update multipliers if convergence criteria are not met
            lambda = lambda - mu * c_val;
            eta = eta / (mu^0.9);
            omega = omega / mu;
  
        else
            % Increase penalty parameter for next iteration
            mu = 100 * mu;
            eta = 1 / (mu^0.1);
            omega = 1 / mu;
        end
    end
    
    % Display final constraint violations
    constraint_violations = calculateC(x, a, D, gamma);
    fprintf('Final constraint violations: ');
    disp(constraint_violations');
end


function cx = calculateC(x, a, D, gamma)
    % Compute the value of the constraints at x
    m = size(a, 1);
    cx = zeros(m, 1);
    for i = 1:m
        cx(i) = a(i,:)*x + gamma*x'*D{i}*x;
    end
end

function obj = computeObjective(x, c, Q0, D0, gamma)
    % Compute the value of the objective at x
    obj = c' * x + x' * Q0 * x + gamma * x' * D0 * x;
end


