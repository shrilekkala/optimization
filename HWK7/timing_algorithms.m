% Script to compare mean computational times to solution for simplex and
% interior point methods, for random linear programs with n = 2m

% Initialize variables
interior_point_meanTimes = [];
interior_point_stdDevTimes = [];
simplex_meanTimes = [];
simplex_stdDevTimes = [];

% Maximum number of iterations allowed for both methods
p = 5000;

% Range of m values to perform the experiments on
m_values = 5:5:250;

% Loop over different values of m
for m = m_values

    fprintf('m = %d \n', m);
    n = 2 * m;

    % Initial seed for reproducibility
    s = 1;      

    % To store times for each instance
    simplex_times = [];
    interior_point_times = [];

    % To check for degeneracy in simplex
    degenerateFlag = false; 
    degeneracyCount = 0;
    
    for instance = 1:10

        % Create a random linear program, 
        % changing random seed for each instance
        [A, b, c] = RandomLinearProgram(m, n, s + instance);

        tic;
        % Solve using the interior point method
        [x_1, obj_1] = longStepInteriorPoint(A, b, c, p);
        interior_point_elapsedTime = toc;
        
        tic;
        % Solve using the simplex method
        [x, obj, isDegenerate, ~, ~] = simplexMethod(A, b, c, p);
        simplex_elapsedTime = toc;
        
        % Check for degeneracy in simplex method
        if isDegenerate
            degeneracyCount = degeneracyCount + 1;
            continue; % Skip this instance
        else
            simplex_times = [simplex_times, simplex_elapsedTime];
        end

        interior_point_times = [interior_point_times, interior_point_elapsedTime];
    end
    
    % Calculate mean and standard deviation of the times
    interior_point_meanTime = mean(interior_point_times);
    interior_point_stdDevTime = std(interior_point_times);
    simplex_meanTime = mean(simplex_times);
    simplex_stdDevTime = std(simplex_times);
    
    % Store results
    interior_point_meanTimes = [interior_point_meanTimes, interior_point_meanTime];
    interior_point_stdDevTimes = [interior_point_stdDevTimes, interior_point_stdDevTime];
    simplex_meanTimes = [simplex_meanTimes, simplex_meanTime];
    simplex_stdDevTimes = [simplex_stdDevTimes, simplex_stdDevTime];
end


% Plot Mean Computation Times
figure;
plot(m_values, simplex_meanTimes, 'b-o', 'LineWidth', 2); hold on;
plot(m_values, interior_point_meanTimes, 'r-s', 'LineWidth', 2);
hold off;
title('Mean Computation Times for Different Problem Sizes');
legend('Simplex Method', 'Interior Point Method');
xlabel('Problem Size (m)');
ylabel('Mean Computation Time (s)');
grid on;

% Plot Standard Deviations
figure; % Create a new figure for the standard deviations
plot(m_values, simplex_stdDevTimes, 'b-o', 'LineWidth', 2); hold on;
plot(m_values, interior_point_stdDevTimes, 'r-s', 'LineWidth', 2);
hold off;
title('Standard Deviation of Computation Times for Different Problem Sizes');
legend('Simplex Method', 'Interior Point Method');
xlabel('Problem Size (m)');
ylabel('Standard Deviation (s)');
grid on;
