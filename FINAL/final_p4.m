% Define the number of instances
m_values = [5:5:50];
soscResults = [];

% Seed for reproducibility
rng(123);

for i = 1:length(m_values)
    fprintf("m = %d \n", m_values(i))

    % Set problem dimensions and parameters
    m = m_values(i);
    n = 2 * m;
    gamma = 10;
    maxIter = 50;

    for j = 1:10
        % Generate a random problem instance with seed j for variability
        [c, Q0, D0, D, a] = RandomNonLinearProgram(m, n, j);
    
        % Solve the problem
        [x_star, lambda, mu] = BCLMethod(c, Q0, D0, D, a, gamma, maxIter);
        
        % Check for SOSC
        soscResults = [soscResults, checkSOSC(x_star, lambda, Q0, D0, D, a, gamma)];
    end
end

% Report frequency
frequencySOSC = sum(soscResults) / length(soscResults);
fprintf('The SOSC is satisfied in %.2f%% of the instances.\n', frequencySOSC * 100);

