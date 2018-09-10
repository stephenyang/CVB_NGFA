
pred = 3;

covZ = eye(K);
for obs = setdiff([1:numgroup], pred)
    covZ = covZ + E_ZWZW{obs};
end

[eigvec, eigval] = eig((covZ+covZ')./2);
invCovZ = eigvec * diag(1./diag(eigval)) * eigvec';

Xpred = zeros(size(Xtest{pred}));
for obs = setdiff([1:numgroup], pred)
    Xpred = Xpred + mean(prec_e{obs}) .* (ZW{pred}' * invCovZ * ZW{obs} * Xtest{obs}')';
end

%%
Xrecon = V * e_ZW{pred};
% ss = (Xrecon - X{pred}).^2;
train_error = zeros(1, size(Xtest{pred}, 2));
for in = 1:size(Xtest{pred}, 2)
    train_error(in) = mean((X{pred}(:, in) - Xrecon(:, in)).^2)/mean(X{pred}(:, in).^2);
end
training_error = mean(train_error);
fprintf('MSE of training data: %.3f.\n', training_error);

%%
error = zeros(1, size(Xtest{pred}, 2));
for in = 1:size(Xtest{pred}, 2)
    error(in) = mean((Xtest{pred}(:, in) - Xpred(:, in)).^2)/mean(Xtest{pred}(:, in).^2);
end
test_error = mean(error);
fprintf('MSE of test data: %.3f.\n', test_error);