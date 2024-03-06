function H = computeHessian(x, lambda, mu, c, Q0, D0, D, a, gamma)
    % Computes the Hessian of the Augmented Lagrangian at a point x
    m = size(a, 1);

    % Hessian of the objective function
    H = 2 * Q0 + 2 * gamma * D0; 

    for i = 1:m
        % Constraint function c_i(x)
        c_i = a(i, :)*x + gamma*x'*D{i}*x;
        % Gradient of the c_i(x)
        grad_c_i = a(i, :)' + 2 * gamma * D{i} * x;
        % Hessian of c_i(x)
        H_c_i = 2 * gamma * D{i};

        % Subtract the Hessian contributions from the equality constraints
        H = H - lambda(i) * H_c_i;

        % Add the Hessian contributions from the quadratic penalty terms
        H = H + mu * (grad_c_i * grad_c_i' + c_i * H_c_i);
    end
end