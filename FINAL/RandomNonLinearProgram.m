function [c, Q0, D0, D, a] = RandomNonLinearProgram(m, n, s)
    % Set the random seed
    rng(s);
    
    % Generate vector c, n x 1, with uniform distribution U[0, 1]
    c = 10 * rand(n, 1);
    
    % Generate vector q0, n x 1, with uniform distribution 0.5 + U[0, 1]
    q0 = 0.5 + rand(n, 1);

    % Create diagonal matrix D0
    D0 = diag(2*rand(n, 1) - 1);
    
    % Create the diagonal matrix Q0 from vector q0
    Q0 = diag(q0);
    
    % Initialize D, a cell array of m diagonal matrices
    D = cell(m, 1);
    
    % Generate m diagonal matrices D_i with elements drawn from U[-1, 1]
    for i = 1:m
        D{i} = diag(2*rand(n, 1) - 1);
    end
    
    % Generate matrix A, m x n, with elements drawn from U[0, 1]
    a = rand(m, n);
    
end
