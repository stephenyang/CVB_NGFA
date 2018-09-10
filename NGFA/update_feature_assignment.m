
%%
Act = zeros(numgroup, K);

%%
tE_n_jk = zeros(numgroup, K);
tV_n_jk = zeros(numgroup, K);
tZ_n_jk = zeros(numgroup, K);

tE_nj_minus_n_jk = zeros(numgroup, K);
tV_nj_minus_n_jk = zeros(numgroup, K);
tZ_nj_minus_n_jk = zeros(numgroup, K);

%%
rand_ig = randperm(numgroup);
for ig = 1:numgroup
    % ig = rand_ig(ig_id);
    Vphi{ig} = V .* repmat((uvec(1, :)' .* prec_e{ig}), 1, K);
    VTV{ig} = Vphi{ig}' * V;
    EVVphi{ig} = sum(EVV .* repmat((uvec(1, :)'.*prec_e{ig}), 1, K));

    %%
%     if iter > 1
%     1./prec_w{ig};
%     epsilon = abs(trace(cov(X{ig}))/Nj(ig) - 1/prec_e{ig}(1)); % + diag(ones(1, Nj(ig))./prec_e{ig}(1))
%     end
    %%
    % rand_ii = randperm(Nj(ig));
    
    for i = 1:Nj(ig)
        % i = rand_ii(ii);
        for k = 1:K
            %%
            E_q_z_ik_eq_one = q_z_ik_eq_one{ig}(i, k);
            V_q_z_ik_eq_one = E_q_z_ik_eq_one * (1 - E_q_z_ik_eq_one);
            
            term1 = G_alpha(ig) * G_beta_k(k) + E_n_jk(ig, k) - E_q_z_ik_eq_one + eps;
            term2 = 0.5 * (V_n_jk(ig, k) - V_q_z_ik_eq_one)/term1/term1;
            term3 = 0.5*( (W{ig}(k, i)^2 + DiagWvar{ig}(k, i)) * EVVphi{ig}(k) - 2 * ( X{ig}(:, i)' * Vphi{ig}(:, k) - ZW{ig}(:, i)' * VTV{ig}(:, k) + VTV{ig}(k, k) * ZW{ig}(k, i) ) * W{ig}(k, i) );
            
            log_q_z_ik_eq_one = log(term1)  -term3 - term2;
            
            %%
            E_q_z_ik_eq_zero = 1 - q_z_ik_eq_one{ig}(i, k);
            V_q_z_ik_eq_zero = E_q_z_ik_eq_zero * (1 - E_q_z_ik_eq_zero);
            
            term4 = G_alpha(ig) * G_one_minus_beta_k(k) + E_nj_minus_n_jk(ig, k) - E_q_z_ik_eq_zero + eps;
            term5 = 0.5 * (V_nj_minus_n_jk(ig, k) - V_q_z_ik_eq_zero)/term4/term4;
            log_q_z_ik_eq_zero = log(term4) - term5;
            
            q_z_ik_eq_one{ig}(i, k) = 1./(1 + exp(log_q_z_ik_eq_zero - log_q_z_ik_eq_one));
               
        end
        
        for k = 1:K
            %%
            E_q_z_idk_eq_one = q_z_ik_eq_one{ig}(i, k);
            V_q_z_idk_eq_one = E_q_z_idk_eq_one * (1 - E_q_z_idk_eq_one);
            Z_q_z_idk_eq_one = log(1- E_q_z_idk_eq_one);
            
            tE_n_jk(ig, k) = tE_n_jk(ig, k) + E_q_z_idk_eq_one;
            tV_n_jk(ig, k) = tV_n_jk(ig, k) + V_q_z_idk_eq_one;
            tZ_n_jk(ig, k) = tZ_n_jk(ig, k) + Z_q_z_idk_eq_one;
            
            %%
            E_q_z_idk_eq_zero = 1- q_z_ik_eq_one{ig}(i, k);
            V_q_z_idk_eq_zero = E_q_z_idk_eq_zero * (1 - E_q_z_idk_eq_zero);
            Z_q_z_idk_eq_zero = log(1- E_q_z_idk_eq_zero + eps);
            
            tE_nj_minus_n_jk(ig, k) = tE_nj_minus_n_jk(ig, k) + E_q_z_idk_eq_zero;
            tV_nj_minus_n_jk(ig, k) = tV_nj_minus_n_jk(ig, k) + V_q_z_idk_eq_zero;
            tZ_nj_minus_n_jk(ig, k) = tZ_nj_minus_n_jk(ig, k) + Z_q_z_idk_eq_zero;
            
        end
                
    end
    
    Z{ig} = binornd(1, q_z_ik_eq_one{ig}');
    
    %%
    elbo_z(iter) = elbo_z(iter) - sum(sum(q_z_ik_eq_one{ig} .* log(q_z_ik_eq_one{ig} + realmin)));
    elbo_z(iter) = elbo_z(iter) + gammaln(G_alpha(ig)) - gammaln(G_alpha(ig) + Nj(ig));
    for k = 1:K
        elbo_z(iter) = elbo_z(iter) ...
            + ((G_alpha_x_beta_k(ig, k) + E_n_jk(ig, k)) > 1e-3) * gammaln(eps + G_alpha_x_beta_k(ig, k) + E_n_jk(ig, k) ) ...
            - (G_alpha_x_beta_k(ig, k) >1e-3 ) * gammaln(eps + G_alpha_x_beta_k(ig, k)) ...
            + ((G_alpha_x_one_minus_beta_k(ig, k) + E_nj_minus_n_jk(ig, k))>1e-3) * gammaln(eps + G_alpha_x_one_minus_beta_k(ig, k) + E_nj_minus_n_jk(ig, k)) ...
            - (G_alpha_x_one_minus_beta_k(ig, k)>1e-3) * gammaln(eps + G_alpha_x_one_minus_beta_k(ig, k));
    end
    %%
end

%%
E_n_jk = tE_n_jk;
V_n_jk = tV_n_jk;
Z_n_jk = tZ_n_jk;

E_nj_minus_n_jk = tE_nj_minus_n_jk;
V_nj_minus_n_jk = tV_nj_minus_n_jk;
Z_nj_minus_n_jk = tZ_nj_minus_n_jk;

%% 
% for ig = 1:numgroup
%     Z{ig} = binornd(1, q_z_ik_eq_one{ig}');
% end


%%
% for ig = 1:numgroup
%     elbo_z(iter) = elbo_z(iter) + gammaln(G_alpha(ig)) - gammaln(G_alpha(ig) + Nj(ig));
%     for k = 1:K
%            elbo_z(iter) = elbo_z(iter) + gammaln( eps + G_alpha_x_beta_k(ig, k) + E_n_jk(ig, k) ) ...
%                - gammaln(eps + G_alpha_x_beta_k(ig, k)) ...
%                + gammaln( eps + G_alpha_x_one_minus_beta_k(ig, k) + E_nj_minus_n_jk(ig, k)) ...
%                - gammaln(eps + G_alpha_x_one_minus_beta_k(ig, k));
%     end
%     
%     elbo_z(iter) = elbo_z(iter) - sum(sum(q_z_ik_eq_one{ig} .* log(q_z_ik_eq_one{ig} + realmin)));
%     
% end
%%