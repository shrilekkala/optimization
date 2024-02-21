% initialize variables
meanTimes = [];
stdDevTimes = [];
degeneracyCounts = [];

% Maximum number of iterations allowed
p = 5000;

% Range of m values to perform the experiments on
m_values = 5:5:250;

% Loop over different values of m
for m = m_values
    fprintf('m = %d \n', m);

    n = 2 * m;

    s = 1; % Initial seed for reproducibility
    times = []; % To store times for each instance
    degenerateFlag = false; % To check for degeneracy
    degeneracyCount = 0;
    
    for instance = 1:10
        % Create a random linear program, 
        % changing random seed for each instance
        [A, b, c] = RandomLinearProgram(m, n, s + instance);

        tic;
        
        % Solve using the simplex method
        [x, obj, isDegenerate, ~, ~] = simplexMethod(A, b, c, p);
        
        elapsedTime = toc;
        
        % Check for degeneracy
        if isDegenerate
            degeneracyCount = degeneracyCount + 1;
            continue; % Skip this instance
        else
            times = [times, elapsedTime];
        end
    end
    
    % Calculate mean and standard deviation of the times
    meanTime = mean(times);
    stdDevTime = std(times);
    
    % Store results
    meanTimes = [meanTimes, meanTime];
    stdDevTimes = [stdDevTimes, stdDevTime];
    degeneracyCounts = [degeneracyCounts, degeneracyCount];
end

% Display results
disp('Mean Times:');
disp(meanTimes);

disp('Standard Deviations:');
disp(stdDevTimes);

disp('Degeneracy Counts:');
disp(degeneracyCounts);

% Plot Mean Computation Times
figure;
plot(m_values, meanTimes, '-o', 'LineWidth', 2);
title('Mean Computation Times for Different Problem Sizes');
xlabel('Problem Size (m)');
ylabel('Mean Computation Time (s)');
grid on;

% Plot Standard Deviations
figure; % Create a new figure for the standard deviations
plot(m_values, stdDevTimes, '-s', 'LineWidth', 2, 'Color', 'r');
title('Standard Deviation of Computation Times for Different Problem Sizes');
xlabel('Problem Size (m)');
ylabel('Standard Deviation (s)');
grid on;
