function alpha_k = backtracking_line_search(fgH, x, pk, alpha_bar, rho, c)
    % Algorithm 3.1
    alpha = alpha_bar;

    [f_next, ~, ~] = fgH(x + alpha * pk);
    [f_current, grad_current, ~] = fgH(x);

    while f_next > f_current + c * alpha * grad_current' * pk
        alpha = rho * alpha; % Contract the step size
        % update f_next using the new alpha:
        [f_next, ~, ~] = fgH(x + alpha * pk);
    end
    alpha_k = alpha;
end
