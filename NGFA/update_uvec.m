
EVnTVn = sum(EVV, 2);

for ig = 1:numgroup
    
    %%
    predX{ig} = V * ZW{ig};
    midval = sum(X{ig} .* (X{ig} - 2 * predX{ig}), 2);
    tmp1 = ZWZW{ig} - diag(diag(ZWZW{ig})) + diag(magWZ{ig});
    %%
    for ni = 1:p
        
        %%
        tmp2 = V(ni, :)' * V(ni, :);
        tmp2 = tmp2 - diag(diag(tmp2)) + diag(EVV(ni, :));
        %%
        alpha_u(1, ni) = 0.5*nu + 0.5*(sum(Nj)) + 0.5*K;
        beta_u(1, ni) = 0.5*nu + 0.5*prec_e{ig}(ni)*(trace(tmp1 * tmp2) + midval(ni)) + 0.5*EVnTVn(ni);
        uvec(1, ni) = alpha_u(1, ni)/beta_u(1, ni);
        
        %%
        if exist('iter')
        elbo_uvec(iter) = elbo_uvec(iter) + 0.5*nu*log(0.5*nu) - gammaln(0.5*nu) ...
             - alpha_u(1, ni) * log(beta_u(1, ni)) + gammaln(alpha_u(1, ni)) ...
             - (0.5*(sum(Nj)) + 0.5*K) * (psi(alpha_u(1, ni)) - log(beta_u(1, ni))) ...
             + (0.5*prec_e{ig}(ni)*(trace(tmp1 * tmp2) + midval(ni)) + 0.5*EVnTVn(ni)) * uvec(1, ni);
        end
    end    
end