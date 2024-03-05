function x_c = computeCauchyPoint(x_k, g_k, B_k, l, u)
    % computes the Cauchy point as per the method outlined in §16.7

    % INPUTS:
    % x_k - Current point in the iteration
    % g_k - Gradient of the objective function at x_k
    % B_k - Approximate Hessian (or actual Hessian) at x_k
    % l - Lower bounds on the variables
    % u - Upper bounds on the variables

   
    n = length(x_k);        % Number of variables
    t_i = zeros(n, 1);      % Initialize breakpoint times
    
    % Calculate breakpoints where each variable hits its bound
    for i = 1:n
        if g_k(i) < 0 && u(i) < Inf
            t_i(i) = (x_k(i) - u(i)) / g_k(i);
        elseif g_k(i) > 0 && l(i) > -Inf
            t_i(i) = (x_k(i) - l(i)) / g_k(i);
        else
            t_i(i) = Inf;  % No breakpoint for inactive variables
        end
    end

    % Sort breakpoints in ascending order and remove duplicates
    t_sorted = unique(sort(t_i));
    
    % Initialize the Cauchy point as the current iterate
    x_c = x_k;
    
    % Iterate over breakpoints to determine the Cauchy point
    for i = 1:length(t_sorted)
        if t_sorted(i) > 0  

            % Calculate the point at the current breakpoint
            x_t = x_k - t_sorted(i) * g_k;

            % Ensure x_t stays within bounds
            x_t = max(min(x_t, u), l);

            % Check for improvement
            if q_k(x_t, x_k, g_k, B_k) >= q_k(x_c, x_k, g_k, B_k)
                break;
            else
                x_c = x_t;  % Update the Cauchy point
            end
        end
    end
end

function val = q_k(x, x_k, g_k, B_k)
    % Evaluate the quadratic model q_k at the point x
    val = (x - x_k)' * g_k + 0.5 * (x - x_k)' * B_k * (x - x_k);
end

