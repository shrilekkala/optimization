function H = rosenbrocknHessian(x)
% function that returns the Hessian of a Rosenbrock function of order n
    n = length(x);

    % Intializations
    H = zeros(n, n);

    for i = 1:(n-1)

        % Hessian calculation
        H(i, i) = H(i, i) + 1200 * x(i)^2 - 400 * x(i+1) + 2;
        if i < n
            H(i, i+1) = H(i, i+1) - 400 * x(i);
            H(i+1, i) = H(i+1, i) - 400 * x(i);
            H(i+1, i+1) = H(i+1, i+1) + 200;
        end
    end
end
