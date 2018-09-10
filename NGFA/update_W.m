
for ig = 1:numgroup
    logdetWsig{ig}=0;
    
    for i = 1:Nj(ig)
        
        %%
        
        B = VTV{ig} - diag(diag(VTV{ig})) + diag(EVVphi{ig}) + eye(K) .* diag((prec_w{ig}(i, :)));
        
        B = diag(q_z_ik_eq_one{ig}(i,:)') * (VTV{ig} - diag(diag(VTV{ig})) + diag(EVVphi{ig})) * diag(q_z_ik_eq_one{ig}(i,:)') + eye(K) .* diag((prec_w{ig}(i, :)));
        
        % B = diag(Z{ig}(:, i)) * (VTV{ig} - diag(diag(VTV{ig})) + diag(EVVphi{ig})) * diag(Z{ig}(:, i)) + eye(K) .* diag((prec_w{ig}(i, :)));

        % [eigvec, eigval] = eig((B + B')./2);
        % midval = eigvec * diag(1./(diag(eigval))) * eigvec';
        midval = pinv(B);
        A = Vphi{ig}' * X{ig}(:, i);
        
        A = diag(q_z_ik_eq_one{ig}(i,:)') * Vphi{ig}' * X{ig}(:, i);
        
        % A = diag(Z{ig}(:, i)) * Vphi{ig}' * X{ig}(:, i);
        
        W_sigma{ig}{i} = midval;
        
        DiagWvar{ig}(:, i) = diag(W_sigma{ig}{i});
        
        W{ig}(:, i) = W_sigma{ig}{i} * A;
        
        EWTW{ig}(i) = W{ig}(:, i)' * W{ig}(:, i) + sum(DiagWvar{ig}(:, i));
        
        logdetWsig{ig}=logdetWsig{ig}+0.5*log(2*pi*(det(W_sigma{ig}{i})+1e-10));
        
    end
    
    ZW{ig} = q_z_ik_eq_one{ig}' .* W{ig}; % Z{ig} .*
    ZWZW{ig} = ZW{ig} * ZW{ig}';
    E_ZWZW{ig} = ZWZW{ig} + W_sigma{ig}{1} .* Nj(ig);
    
    %%
    if strcmp(sparsity_w, 'factor-wise')
        e = e0 + 0.5 * Nj(ig) * ones(Nj(ig), K);
        % f = f0 + (Z{ig}.*W{ig}).*(Z{ig}.*W{ig}) + (Z{ig}.*DiagWvar{ig});
        f = f0 + repmat(0.5 * (diag(W{ig} * W{ig}') + sum(DiagWvar{ig}, 2))', Nj(ig), 1);
        % f = f0 + (W{ig}).*(W{ig}) + (DiagWvar{ig});
    elseif strcmp(sparsity_w, 'element-wise')
        e = e0 + 0.5;
        f = f0 + ( (q_z_ik_eq_one{ig}'.*W{ig}).*(q_z_ik_eq_one{ig}'.*W{ig}) + (q_z_ik_eq_one{ig}'.*DiagWvar{ig}) )';
    end
    prec_w{ig} = e./f;
    
    
    %%
    if strcmp(sparsity_w, 'factor-wise')
        %%
        elbo_w(iter) = elbo_w(iter) + 0.5* Nj(ig)*sum(log(prec_w{ig}(1, :))) - 0.5*trace( diag(prec_w{ig}(1, :)) .* diag((diag(W{ig} * W{ig}') + sum(DiagWvar{ig}, 2))) ) ...
            + logdetWsig{ig};
        %%
        ee = e0 + 0.5 * Nj(ig);
        ff = f0 + 0.5 * (diag(W{ig} * W{ig}') + sum(DiagWvar{ig}, 2));
        elbo_prec_w(iter) = elbo_prec_w(iter) + sum(sum(e0*log(f0) - gammaln(e0) - ee.*log(ff) + gammaln(ee)+ (e0-ee)*(psi(ee) - log(ff)) + (1-f0./ff)*ee) );
    elseif strcmp(sparsity_w, 'element-wise')
        elbo_w(iter) = elbo_w(iter) + 0.5 * sum(sum( psi(e) - log(f) )) - 0.5*sum(sum(prec_w{ig} .* (W{ig}.* W{ig} + DiagWvar{ig})' ) )...
            +logdetWsig{ig};
        
        elbo_prec_w(iter) = elbo_prec_w(iter) + sum(sum(e0*log(f0) - gammaln(e0) - e.*log(f) + gammaln(e)+ (e0-e).*(psi(e) - log(f)) + (1-f0./f).*e) );
    end
    %%
    
end