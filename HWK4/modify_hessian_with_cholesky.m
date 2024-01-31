function Bk = modify_hessian_with_cholesky(H, beta)
    % Algorithm 3.3
    if min(diag(H)) > 0
        tau = 0;
    else
        tau = -min(diag(H)) + beta;
    end
    
    k = 0;
    while true
        % apply Cholesky algorithm
        [L, p] = chol(H + tau * eye(size(H)), 'lower');
        % p = 0 signifies that the factorization was successful
        if p == 0
            Bk = L * L';
            return;
        else
            tau = max(2 * tau, beta);
        end
        k = k + 1;
        % exit condition
        if k > 100
            error('Hessian could not be modified to be positive definite.');
        end
    end
end