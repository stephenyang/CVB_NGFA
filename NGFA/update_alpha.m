  
E_log_u_j = zeros(numgroup, K);

E_m_jk = zeros(numgroup, K);

E_s_jk = zeros(numgroup, K);

for ig = 1:numgroup
    for k = 1:K
        %%
        if  E_n_jk(ig, k) > 1
            Pplus_n_jk = 1 - exp(Z_n_jk(ig, k));
            Eplus_n_jk = E_n_jk(ig, k)/Pplus_n_jk;
            Vplus_n_jk = V_n_jk(ig, k)/Pplus_n_jk - exp(Z_n_jk(ig, k)) * Eplus_n_jk;
            psi1 = psi(G_alpha_x_beta_k(ig, k) + Eplus_n_jk);
            psi3 = psi(2, G_alpha_x_beta_k(ig, k) + Eplus_n_jk);
            % psi3 = trigamma(G_alpha_x_beta_k(k) + Eplus_n_jk);
            E_m_jk(ig, k) = G_alpha_x_beta_k(ig, k) * Pplus_n_jk * (psi1 - digamma_G_alpha_x_beta_k(ig, k) + Vplus_n_jk * psi3/2);
        else
            E_m_jk(ig, k) =0;
        end
        
        %%
        if  (Nj(ig) - E_n_jk(ig, k)) > 1
            Pplus_nj_minus_n_jk = 1 - exp(Z_nj_minus_n_jk(ig, k));
            Eplus_nj_minus_n_jk = E_nj_minus_n_jk(ig, k)/Pplus_nj_minus_n_jk;
            Vplus_nj_minus_n_jk = V_nj_minus_n_jk(ig, k)/Pplus_nj_minus_n_jk - exp(Z_nj_minus_n_jk(ig, k)) * Eplus_nj_minus_n_jk;
            psi_a = psi(G_alpha_x_one_minus_beta_k(ig, k) + Eplus_nj_minus_n_jk);
            psi_b = psi(2, G_alpha_x_one_minus_beta_k(ig, k) + Eplus_nj_minus_n_jk);
            % psi_b = trigamma(G_alpha_x_one_minus_beta_k(k) + Eplus_nj_minus_n_jk);
            E_s_jk(ig, k) = G_alpha_x_one_minus_beta_k(ig, k) * Pplus_nj_minus_n_jk ...
                * (psi_a - digamma_G_alpha_x_one_minus_beta_k(ig, k) + Vplus_nj_minus_n_jk * psi_b/2);
        else
            E_s_jk(ig, k) = 0;
        end
    end
    E_log_u_j(ig, 1:K) = psi(E_alpha(ig) +1) - psi(E_alpha(ig) + Nj(ig)+1);
end

E_m = sum(E_m_jk, 2);
E_s = sum(E_s_jk, 2);

E_l_jk = zeros(numgroup, K);
for ig = 1:numgroup
    E_l_jk(ig, :) = Nj(ig)/(Nj(ig) + G_alpha(ig));
end
E_l = sum(E_l_jk, 2);

E_log_u = sum(E_log_u_j, 2);

for ig = 1:numgroup
aa(ig) = a_alpha(ig) + E_m(ig) + E_s(ig) - E_l(ig);
bb(ig) = b_alpha(ig) - E_log_u(ig);

E_alpha(ig) = aa(ig)/bb(ig);
G_alpha(ig) = exp(psi(aa(ig) + eps))/bb(ig);

elbo_alpha(iter) = elbo_alpha(iter) + a_alpha(ig) * log(b_alpha(ig)) - gammaln(a_alpha(ig)) - aa(ig)*log(bb(ig)) + gammaln(aa(ig))...
    +(bb(ig) - b_alpha(ig)) * G_alpha(ig) + (a_alpha(ig) - aa(ig)) * (psi(aa(ig)) - log(bb(ig)));
end

%%
% elbo_alpha(iter) = a_alpha * log(b_alpha) - gammaln(a_alpha) - a*log(b) + gammaln(a)...
%     +(b - b_alpha) * G_alpha + (a_alpha - a) * (psi(a) - log(b));
% for ig = 1:numgroup
% elbo_alpha(iter) = elbo_alpha(iter) + a_alpha(ig) * log(b_alpha(ig)) - gammaln(a_alpha(ig)) - a(ig)*log(b(ig)) + gammaln(a(ig))...
%     +(b(ig) - b_alpha(ig)) * E_alpha(ig) + (a_alpha(ig) - a(ig)) * (psi(a(ig)) - log(b(ig)));
% end
%%

