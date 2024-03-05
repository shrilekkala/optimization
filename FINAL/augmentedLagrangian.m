function LA = augmentedLagrangian(x, lambda, mu, c, Q0, D0, D, a, gamma)
    % Computes the value of the Augmented Lagrangian at a point x

    f_x = c'*x + x'*Q0*x + gamma*x'*D0*x;

    constraint_penalty = 0;

    for i = 1:length(lambda)
        ci_x = a(i,:)*x + gamma*x'*D{i}*x;
        constraint_penalty = constraint_penalty - lambda(i)*ci_x + (mu/2)*(ci_x)^2;
    end
    
    LA = f_x + constraint_penalty;
end