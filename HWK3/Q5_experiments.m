% initializations
max_n = 1000;
times = zeros(max_n, 2);
tol = 1e-6; % stopping criterion based on gradient norm
max_iter = 100; % maximum iterations
initial_points = [5, 10]; % initial points to test

for i = [1,2]
    initial_point = initial_points(i)
    fprintf('Running experiments for x_0 with all entries equal to %d\n', initial_point);
    for n = 1:max_n
        x0 = initial_point * ones(n, 1);

        tic; % start timer
        x_star = NMHM(@rosenbrocknfgH, x0, tol, max_iter);
        times(n, i) = toc; % stop timer and store compute time

        % Check if total compute time for the run exceeds 2 minutes
        if sum(times(:, i)) > 120
            fprintf('Compute time exceeded 1 minute at n = %d\n', n);
            break;
        end

        % Output the results
        fprintf('n = %d, Time = %f seconds\n', n, times(n));
    end
end

% Plot the compute times for the entry of all 10s
figure;
plot(1:n, times(1:n, 2), 'LineWidth', 2);
title(sprintf('Compute times for x_0 with all entries 10'));
xlabel('Order n of Rosenbrock function');
ylabel('Compute time (s)');
grid on;

