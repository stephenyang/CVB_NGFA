
%% prerequisites
% verbose
% K
% isotropic
% modeltype
% robust_noise
% robust_factor

%% initialization

get_default_setting;

initialization;

if strcmp(verbose, 'true')
    figure(123);
    for ig = 1:numgroup
        subplot(numgroup+1, 3, 1+(ig-1)*3);
        imagesc(AA{ig});colormap('gray');axis off;
        title(['True factor loadings', num2str(ig)])
    end
end

%% main loop
tic;
for iter = 1:maxiter
    
    %%
    update_feature_assignment;
    
    update_alpha;
    
    update_beta;
    
    %%

    update_W;

    
    update_V;
    
    for ig = 1:numgroup
        % e_Z{ig} = binornd(1, Z{ig});
        e_Z{ig} =Z{ig};
        e_ZW{ig} = e_Z{ig} .* W{ig};
        % e_Z{ig} =Z{ig};
    end
    
    update_prec_e;
    
    if strcmp(modeltype, 'robust')
        update_uvec;
        
        update_nv;
    end
    
    %% elbo
    elbo(iter) = elbo_x(iter)  + elbo_w(iter) + elbo_z(iter) + elbo_v(iter) ...
    + elbo_prec_w(iter) + elbo_prec_e(iter) ...
    + elbo_alpha(iter) + elbo_beta(iter) + elbo_uvec(iter) ;
    
    %%
    zsum = sum(e_Z{1}, 2);
    for ig = 2:numgroup
        zsum = zsum + sum(e_Z{ig}, 2);
    end
    active_ind = find(zsum>0);
    Kactive = numel(find(zsum>0));
    %%
    if strcmp(verbose, 'true')
        figure(123);
        [~, ord] = sort(zsum,  'descend');
        for ig = 1:numgroup
            subplot(numgroup+1, 3, 2+(ig-1)*3);
            imagesc(e_ZW{ig}(active_ind, :)');colormap('gray');title(['Inferred factor loadings', num2str(ig)]);axis off;
            subplot(numgroup+1, 3, 3+(ig-1)*3)
            imagesc(e_Z{ig}([active_ind' setdiff(ord, active_ind)'], :)');colormap('gray');title(['Inferred binary matrix', num2str(ig)]);axis off;
            drawnow;
        end
        %%
        if iter > 1
            subplot(numgroup+1, 3, ((numgroup+1)*3-2):((numgroup+1)*3));plot([iter-1 iter],[elbo(iter-1) elbo(iter)],'b-');hold on;xlim([-5 miniter + 5]); %ylim([elbo(1) elbo(iter)+1e4]);
            ylabel('ELBO');xlabel('epoch');
        end
    disp(['iterations:' num2str(iter) ' ELBO: ' num2str(elbo(iter))]);    
    end
    
    if iter > miniter
        % abs((elbo(iter) - elbo(iter-1))./elbo(iter-1)) 
        if abs((elbo(iter) - elbo(iter-1))./elbo(iter-1)) < 1e-3
            return;
        end
    end
    
end
toc