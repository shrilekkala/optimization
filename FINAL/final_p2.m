% Define m values and preallocate arrays for results
m_values = [5, 10:10:100];
mean_times = zeros(length(m_values), 1);
std_times = zeros(length(m_values), 1);

% Seed for reproducibility
rng(2468);

% Loop over each m value
for i = 1:length(m_values)
    fprintf("Solving problems for m = %d \n", i)

    % Set problem dimensions and parameters
    m = m_values(i);
    n = 2 * m;
    gamma = 0;
    maxIter = 100;
    
    % Preallocate an array for storing times
    times = zeros(10, 1);
    
    for j = 1:10
        % Generate a random problem instance with seed with j for variability
        [c, Q0, D0, D, a] = RandomNonLinearProgram(m, n, j);
        
        % Time the BCLMethod for this instance
        tic;
        x_star = BCLMethod(c, Q0, D0, D, a, gamma, maxIter);
        times(j) = toc;
    end
    
    % Compute mean and standard deviation of times
    mean_times(i) = mean(times);
    std_times(i) = std(times);
end

% Plot mean times
figure;
plot(m_values, mean_times, '-o');
title('Mean Times');
xlabel('m value');
ylabel('Mean Time (seconds)');
grid on;

% Plot standard deviation times
figure;
plot(m_values, std_times, '-o');
title('Standard Deviation of Times');
xlabel('m value');
ylabel('Standard Deviation (seconds)');
grid on;
