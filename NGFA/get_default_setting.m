
%%
numgroup = numel(X);
Nj = zeros(1, numgroup);
for ig = 1:numgroup
    Nj(ig) = size(X{ig}, 2);
end
p = size(X{1}, 1);

%%
modeltype = 'robust';
robust_noise = 'true';
robust_factor = 'false';
isotropic = 'true';
K = min([Nj p]);
miniter  = 2 * 1e2;
maxiter = 1e3;
verbose = 'true';
sparsity_w = 'element-wise';

%% pr of nv
nu = 1;

%% pr of prec_v
c0 = 1* ones(p, K); % 10^-6 * ones(p, K);
d0 = 1* ones(p, K); % 10^-6 * ones(p, K);

%% pr of prec_w
e0 = 1;%1e-15;
f0 = 1;%1e-15;

%% pr of prec_e
if ~strcmp(isotropic, 'true')
        g0 = 1;%e-10 * ones(p, 1);
        h0 = 10;%e-10* ones(p, 1);
else
    g0 = 1;%1e-15;
    h0 = 10;% 1e-15* ones(p, 1);
end

%% gamma global concentration parameter
eta = 1;
a_gamma = 1;
b_gamma = eta;
E_gamma_a = a_gamma/b_gamma;
E_gamma_b = eta;