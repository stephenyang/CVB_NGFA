
%%
% modeltype = 'robust';
% robust_noise = 'true';
% robust_factor = 'false';
% isotropic = 'true';
% K = 50;
% miniter = 1e2;
% maxiter = 1e5;
% verbose = 'true';
% sparsity_w = 'element-wise';
% 
% %%
% numgroup = numel(X);
% 
% Nj = zeros(1, numgroup);
% for ig = 1:numgroup
%     Nj(ig) = size(X{ig}, 2);
% end
% 
% p = size(X{1}, 1);

Kactive = K;
active_ind = 1:K;

%%
optimal_mse = 1e3;
elbo_x = zeros(1, maxiter);
elbo_w = zeros(1, maxiter);
elbo_z = zeros(1, maxiter);
elbo_v = zeros(1, maxiter);
elbo_prec_w = zeros(1, maxiter);
elbo_prec_e = zeros(1, maxiter);
gg0 = 1* ones(p, 1);
gg = gg0 + 0.5*Nj(ig);
elbo_alpha = zeros(1, maxiter);
elbo_beta = zeros(1, maxiter);
elbo_uvec = zeros(1, maxiter);
elbo = zeros(1, maxiter);
elbo_trace = zeros(1, maxiter);    

%%
nu = 10;
alpha_u = zeros(1, p);
beta_u = zeros(1, p);
uvec = ones(1, p);

%% 
c0 = 1* ones(p, K); % 10^-6 * ones(p, K); % prior of prec_v
d0 = 1* ones(p, K); % 10^-6 * ones(p, K);

e0 = 1;%10^-6; % prior of prec_w
f0 = 1;%10^-6;

if ~strcmp(isotropic, 'true')
   % if strcmp(noise_structure, 'sample_noise')
        g0 = 1 * ones(p, 1); %  10^-6 * ones(p, 1);
        h0 = 10* ones(p, 1);
   % elseif strcmp(noise_structure, 'variable_noise')
   %     %%
   % end
else
    g0 = 1;
    h0 = 10* ones(p, 1);
end
g = g0; 
h = h0;

prec_v = ones(p, K);


for ig = 1:numgroup
    
    if p > Nj(ig)
        [Xu, Xlambda, Xv] = svd(X{ig}, 'econ');
    else
        [Xu, Xlambda, Xv] = svd(X{ig});
    end
    
    Z{ig} = double(rand(K, Nj(ig)) > (eps + abs(1-p/Nj(ig))));
    
    q_z_ik_eq_one{ig} = Z{ig}';
    
    if K <= Nj(ig)
        W{ig} = Xv(1:K, :);
    else
        W{ig} = [Xv(1:Nj(ig), :);randn(K- Nj(ig), Nj(ig))];
    end
    % W{ig} = [randn(K, Nj(ig))];
    
    ZW{ig} = Z{ig} .* W{ig};
    ZWZW{ig} = ZW{ig} * ZW{ig}';
    
    DiagWvar{ig} = zeros(K, Nj(ig));
    
    if strcmp(sparsity_w, 'factor-wise')
        e = e0 * ones(Nj(ig), K);% + 0.5 * Nj(ig);
        f = f0 * ones(Nj(ig), K);% + 0.5 * (diag(W{ig} * W{ig}') + sum(DiagWvar{ig}, 2));
    elseif strcmp(sparsity_w, 'element-wise')
        e = e0 + 0.5* ones(Nj(ig), K);
        f = f0 + ( (Z{ig}.*W{ig}).*(Z{ig}.*W{ig}) + (Z{ig}.*DiagWvar{ig}) )';
    end
%     e = e0;
%     f = f0;%+(Z{ig}.*W{ig}).*(Z{ig}.*W{ig}) + (Z{ig}.*DiagWvar{ig});
    prec_w{ig} = e./f;

    % prec_w{ig} = e0/f0 * ones(K, Nj(ig));
    prec_e{ig} = g0./h0; %ones(p, 1);

end

V = randn(p, K);
V_sigma = (1e-6) * ones(p, K);
update_V;

%% alpha
a_alpha = 1 * ones(1, numgroup);
b_alpha = 1 * ones(1, numgroup);

E_alpha = a_alpha./b_alpha;
G_alpha = exp(psi(a_alpha))./b_alpha;

%% gamma
eta = 1;
a_gamma = 1;
b_gamma = eta;
E_gamma_a = a_gamma/b_gamma;
E_gamma_b = eta;

%% beta
E_beta_k = E_gamma_a/K * ones(1, K);
E_one_minus_beta_k = (E_gamma_b) * (K-1)/K * ones(1, K);
a = 1/K * ones(1, K);
b = (K-1)/K * ones(1, K);
G_beta_k = E_gamma_a/K * ones(1, K); % exp(psi(a) - psi(a+b));
G_one_minus_beta_k = (E_gamma_b) * (K-1)/K * ones(1, K); % exp(psi(b) - psi(a+b));

%% auxiliary variables
E_log_u = 0;

E_m_jk = zeros(numgroup, K);
E_m_k = zeros(1, K);

E_s_jk = zeros(numgroup, K);
E_s_k = zeros(1, K);

G_alpha_x_beta_k = zeros(numgroup, K);

digamma_G_alpha_x_beta_k = zeros(numgroup, K);

%%
initialize_sufficient_statistics;
    
if strcmp(modeltype, 'robust')
    update_prec_e;
    update_uvec;
    update_nv;
end
