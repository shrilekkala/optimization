function [tau] = find_tau(pU, pB, delta_k)
    % Finds the scalar tau that minimizes the model along the dogleg path within the trust region
    % i.e. finds the solution to equation at the bottom of p75
    a = (pB - pU)' * (pB - pU);
    b = 2 * pU' * (pB - pU);
    c = pU' * pU - delta_k^2;
    
    % Solve the quadratic equation a*x^2 + b*x + c = 0, note here x=tau-1
    x = roots([a, b, c]);
    
    % We are only interested in the real and positive roots within [0, 1]
    % So tau is be between 1 and 2
    x = x(imag(x) == 0 & x >= 0 & x <= 1);
    
    % If there are two such roots, choose the one that gives the minimum m_k(p)
    if length(x) > 1
        disp(6)
        % Choose the tau that gives the smallest norm, which is equivalent to the minimum m_k(p)
        [~, idx] = min(abs(x + 1));
        x = x(idx);
    end

    tau = x + 1;
end