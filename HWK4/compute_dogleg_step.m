function [p_k] = compute_dogleg_step(g_k, B_k, delta_k)
    % Computes the dogleg step p_k for a given gradient g_k,
    % Hessian B_k, and trust region radius delta_k
    
    % Compute the Cauchy point pU
    pU = - (g_k' * g_k) / (g_k' * B_k * g_k) * g_k;
    
    % Check if the full Newton step is within the trust region
    if norm(pU) >= delta_k
        
        % The Cauchy point is outside the trust region, so return the scaled Cauchy point
        p_k = (delta_k / norm(pU)) * pU;
    else
        % Compute the full Newton step (Assuming B_k is positive definite)
        pB = - B_k \ g_k;
        
        if norm(pB) <= delta_k
            % The full Newton step is within the trust region
            p_k = pB;

        else
            % The full Newton step is outside the trust region, so find the intersection
            % Solve the scalar quadratic equation ||pU + tau*(pB - pU)||^2 = delta_k^2
            tau = find_tau(pU, pB, delta_k);
            p_k = pU + (tau - 1) * (pB - pU);
        end
    end
end
