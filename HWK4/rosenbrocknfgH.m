function [f, g, H] = rosenbrocknfgH(x)
    n = length(x);

    % Intializations
    f = 0;
    g = zeros(n, 1);
    H = zeros(n, n);

    for i = 1:(n-1)
        f = f + 100 * (x(i+1) - x(i)^2)^2 + (x(i) - 1)^2;

        % Note each x_i has contributions to gradient and Hessian
        % at most from x_(i-1) and x(i+1) and no other variables
        
        % Gradient calculation
        g(i) = g(i) - 400 * x(i) * (x(i+1) - x(i)^2) + 2 * (x(i) - 1);
        if i < n
            g(i+1) = g(i+1) + 200 * (x(i+1) - x(i)^2);
        end
        
        % Hessian calculation
        H(i, i) = H(i, i) + 1200 * x(i)^2 - 400 * x(i+1) + 2;
        if i < n
            H(i, i+1) = H(i, i+1) - 400 * x(i);
            H(i+1, i) = H(i+1, i) - 400 * x(i);
            H(i+1, i+1) = H(i+1, i+1) + 200;
        end
    end
end
