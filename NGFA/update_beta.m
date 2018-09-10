
E_m_k = sum(E_m_jk, 1);

E_s_k = sum(E_s_jk, 1);

for k = 1:K
    a = E_gamma_a/K + E_m_k(k);
    b = E_gamma_b * (K - 1)/K + E_s_k(k);
    
    psi_a = psi(a);
    psi_b = psi(b);
    psi_a_plus_b = psi(a+b);
    
    E_beta_k(k) = a/(a+b);
    G_beta_k(k) = exp(psi_a - psi_a_plus_b);
    
    E_one_minus_beta_k(k) = b/(a+b);
    G_one_minus_beta_k(k) = exp(psi_b - psi_a_plus_b);
    
    %%
    elbo_beta(iter) = elbo_beta(iter) + gammaln(1) - gammaln(1./K) - gammaln(1-1/K) - gammaln(1 + E_m_k(k) + E_s_k(k)) + gammaln(1/K + E_m_k(k)) + gammaln( (1 - 1/K) + E_s_k(k)) ...
        - E_m_k(k) * (psi_a - psi_a_plus_b) - E_s_k(k) * (psi_b - psi_a_plus_b);
    %%
end

