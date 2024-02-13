function [alpha_k, n_evals] = backtracking_line_search(fgH, x, pk, alpha_bar, rho, c)
    % Count the number of function evaluations
    n_evals = 0;

    % Algorithm 3.1
    alpha = alpha_bar;

    [f_next, ~, ~] = fgH(x + alpha * pk);
    [f_current, grad_current, ~] = fgH(x);

    % Update function evluations count
    n_evals = n_evals + 2;

    while f_next > f_current + c * alpha * grad_current' * pk
        alpha = rho * alpha; % Contract the step size
        % update f_next using the new alpha:
        [f_next, ~, ~] = fgH(x + alpha * pk);

        % Update function evluations count
        n_evals = n_evals + 1;
    end
    alpha_k = alpha;
end
