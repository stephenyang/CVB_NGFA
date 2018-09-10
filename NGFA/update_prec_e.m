

for ig = 1:numgroup
    g = g0 + 0.5 * Nj(ig);
    predX{ig} = V * ZW{ig};
    midval = sum(X{ig} .* (X{ig} - 2 * predX{ig}), 2);
    tmp1 = ZWZW{ig} - diag(diag(ZWZW{ig})) + diag(magWZ{ig});
    for n = 1:p
        tmp2 = V(n, :)' * V(n, :);
        tmp2 = tmp2 - diag(diag(tmp2)) + diag(EVV(n, :));
        h(n) = h0(n) + 0.5 * uvec(1, n) * (trace(tmp1 * tmp2) + midval(n));
    end
    
    if ~strcmp(isotropic, 'true')
        prec_e{ig} = g./h;
    else % if strcmp(noise_structure, 'sample_noise')
        g = g0 + 0.5*Nj(ig) * p;
        sum_h = sum(h);
        prec_e{ig} = g/sum_h * ones(p, 1);
%     elseif strcmp(noise_structure, 'variable_noise')
        %%
        
%     elseif strcmp(noise_structure, 'element_noise')
%         %%
    end
    
    %%
    if exist('iter')
        if strcmp(isotropic, 'true')
            elbo_x(iter) = elbo_x(iter) + 0.5*p*Nj(ig)*log(2*pi) + 0.5*p*Nj(ig)*(psi(g)-log(sum_h))-g/sum_h * (sum_h-sum(h0));
        else % if strcmp(noise_structure, 'sample_noise')
            elbo_x(iter) = elbo_x(iter) + 0.5*p*Nj(ig)*log(2*pi) + sum( 0.5*Nj(ig)*(psi(g)-log(h))-(g./h) .* (h-sum(h0)) );
        end
        elbo_prec_e(iter) = elbo_prec_e(iter) -(gg'*log(h) - gg0'*log(h0) - sum(gammaln(gg)) + sum(gammaln(gg0)) + (gg-gg0)'*(psi(g)-log(h)) - gg'*(1-h0./h) );
    end
    %%    
end

