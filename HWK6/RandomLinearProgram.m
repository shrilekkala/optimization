function [A, b, c] = RandomLinearProgram(m, n, s)
    % Matrix Dimensions: m, n
    % Random seed: s

    % Set random seed
    rng(s);

    % Matrix with Standard Normal entries
    A = randn(m, n);
    
    % Vector with first m components 1 and the rest 0
    e = [ones(m, 1); zeros(n-m, 1)];
    
    % Compute b
    b = A * e;
    
    % Create a random objective vector c
    u = rand(n, 1); % Vector with each entry uniformly distributed in (0,1)
    c = [ones(n, 1) + 50*u];
end
