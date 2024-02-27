function alpha = lineSearchAlgorithm(x, s, delta_x, delta_s, gamma, mu)
    % Implements a backtracking line search algorithm to find the maximum 
    % alpha_k such that the next iterate is in N_{-infinity}(gamma)

    alpha = 1;
    while true
        new_x = x + alpha * delta_x;
        new_s = s + alpha * delta_s;
        if all(new_x > 0) && all(new_s > 0) && all(new_x .* new_s >= gamma * mu)
            break;
        end
        alpha = alpha * 0.9; % Decrease alpha_k
        if alpha < 1e-4      % Threshold to avoid infinite loops
            break;
        end
    end
end