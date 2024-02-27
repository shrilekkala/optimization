function [A, b, c] = RandomLinearProgram(m, n, s)
    % Matrix Dimensions: m, n
    % Random seed: s

    % Set random seed
    rng(s);

    % Matrix with Standard Normal entries
    A = randn(m, n);
    
    % nx1 vector of ones
    e = ones(n, 1);
    
    % Compute b
    b = A * e;
    
    % Create a random objective vector c
    u = rand(n, 1); % nx1 vector with entries uniformly distributed in (0,1)
    c = e + 100*u;
end
