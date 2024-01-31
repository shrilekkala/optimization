function m_k = model_f(f_k, g_k, B_k, p)
    % (4.3)
    % Compute the model function m_k(p) = f_k + g_k^T p + 1/2 p^T B_k p
    m_k = f_k + g_k' * p + 0.5 * p' * B_k * p;
end
