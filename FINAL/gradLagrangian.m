function gradL = gradLagrangian(x, lambda, mu, c, Q0, D0, D, a, gamma)
    % Computes the Hessian of the Augmented Lagrangian at a point x

    % Calculate the gradient of the objective function
    grad_f_x = c + 2 * Q0 * x + 2 * gamma * D0 * x;
    
    % Initialize the gradient of the constraint penalty
    grad_c_penalty = zeros(length(x), 1);
    
    % Calculate the gradient of the constraints weighted by lambda and mu
    for i = 1:length(lambda)
        grad_c_i = a(i, :)' + 2 * gamma * D{i} * x;
        c_i = a(i, :)*x + gamma*x'*D{i}*x;
        grad_c_penalty = grad_c_penalty - lambda(i) * grad_c_i + mu * c_i * grad_c_i;
    end
    
    % The gradient of the augmented Lagrangian
    gradL = grad_f_x + grad_c_penalty;
end