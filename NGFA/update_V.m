
sum_magWZphi = zeros(K, p);
for ig = 1:numgroup
    magWZ{ig} = sum(q_z_ik_eq_one{ig}' .* (W{ig} .^ 2 + DiagWvar{ig}), 2);
    sum_magWZphi = sum_magWZphi + magWZ{ig} * (uvec(1, :)' .* prec_e{ig})';
end

for k = 1:K
    V_sigma(:, k) = 1./(sum_magWZphi(k, :)' + strcmp(robust_factor, 'false') * (prec_v(:, k)));% + strcmp(robust_factor, 'true') * (uvec(1, :)' .* prec_v(:, k)) );
    
    V(:, k) = 0;
    
    sum_phi_Xmid = zeros(p, 1);
    for ig = 1:numgroup
        Xmid{ig} = X{ig} * ZW{ig}(k, :)' - V * ZWZW{ig}(:, k);
        phi_Xmid{ig} = (uvec(1, :)' .* prec_e{ig}) .* Xmid{ig};
        sum_phi_Xmid = sum_phi_Xmid + phi_Xmid{ig};
    end
    
    V(:, k) = V_sigma(:, k) .* sum_phi_Xmid;
 
end

EVV = V.^2 + V_sigma;

%% prec_V
c = c0 + 0.5;
d = d0 + 0.5*EVV;
% prec_v = c./d;

%%
if exist('iter')
elbo_v(iter) = - 0.5*p*K*log(2*pi) + 0.5*sum(sum(2*pi*abs(V_sigma)))...
    -0.5*p*K*log(2*pi) + 0.5*sum(sum(psi(c)-log(d))) - 0.5 * sum(sum((c./d).*(d-d0)));
end
%%
