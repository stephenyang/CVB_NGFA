
%%
E_n_jk = zeros(numgroup, K);
V_n_jk = zeros(numgroup, K);
Z_n_jk = zeros(numgroup, K);

E_nj_minus_n_jk = zeros(numgroup, K);
V_nj_minus_n_jk = zeros(numgroup, K);
Z_nj_minus_n_jk = zeros(numgroup, K);

for ig = 1:numgroup
    Vphi{ig} = V .* repmat(prec_e{ig}, 1, K);
    
    VTV{ig} = Vphi{ig}' * V;
    
    EVVphi{ig} = sum(EVV .* repmat(prec_e{ig}, 1, K));
    
    for i = 1:Nj(ig)
        for k = 1:K
            %%
            E_q_z_idk_eq_one = q_z_ik_eq_one{ig}(i, k);
            V_q_z_idk_eq_one = E_q_z_idk_eq_one * (1 - E_q_z_idk_eq_one);
            Z_q_z_idk_eq_one = log(1- E_q_z_idk_eq_one);
            
            E_n_jk(ig, k) = E_n_jk(ig, k) + E_q_z_idk_eq_one;
            V_n_jk(ig, k) = V_n_jk(ig, k) + V_q_z_idk_eq_one;
            Z_n_jk(ig, k) = Z_n_jk(ig, k) + Z_q_z_idk_eq_one;
            
            %%
            E_q_z_idk_eq_zero = 1- q_z_ik_eq_one{ig}(i, k);
            V_q_z_idk_eq_zero = E_q_z_idk_eq_zero * (1 - E_q_z_idk_eq_zero);
            Z_q_z_idk_eq_zero = log(1- E_q_z_idk_eq_zero);
            
            E_nj_minus_n_jk(ig, k) = E_nj_minus_n_jk(ig, k) + E_q_z_idk_eq_zero;
            V_nj_minus_n_jk(ig, k) = V_nj_minus_n_jk(ig, k) + V_q_z_idk_eq_zero;
            Z_nj_minus_n_jk(ig, k) = Z_nj_minus_n_jk(ig, k) + Z_q_z_idk_eq_zero;
        end
    end
end


for ig = 1:numgroup
    for k = 1:K
        G_alpha_x_beta_k(ig, k) = G_alpha(ig) * G_beta_k(k);
        digamma_G_alpha_x_beta_k(ig, k) = psi(G_alpha_x_beta_k(ig, k) + eps);
        
        G_alpha_x_one_minus_beta_k(ig, k) = G_alpha(ig) * G_one_minus_beta_k(k);
        digamma_G_alpha_x_one_minus_beta_k(ig, k) = psi(G_alpha_x_one_minus_beta_k(ig, k) + eps);
    end
end
