
log_u = psi(alpha_u) - log(beta_u);

delta = 1/sum(Nj) * sum(sum(log_u - uvec));

goalfun = @(x) (1 + log(x/2) - psi(x/2) + delta);

nu = fsolve(goalfun, nu, optimset('Display','off'));

% elbo_trace(iter) = nu;