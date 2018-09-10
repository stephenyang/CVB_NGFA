clear;close all;clc
addpath(genpath('NGFA/'));

for nsamp = 20
for seed = 1
    %%
    file = ['data/BASS_sim3_n',num2str(nsamp),'_seed',num2str(seed),'.mat'];
    
    load(file);% X{2} = X{4};X = {X{1};X{2}};
    numgroup = numel(X);
    % rng(1);
    for rep = 1
        ssi = zeros(1, numgroup);
        %%
        main_cvb_NGFA;
        
        %%
        zsum = sum(Z{1}, 2);
        for ig = 2:numgroup
            zsum = zsum + sum(Z{ig}, 2);
        end
        [~, ord] = sort(zsum,  'descend');
        act_factor_ind = find(zsum>1e-10);
        activeK = numel(act_factor_ind);
        
        %%
        for ig = 1:numgroup
            ssi(ig) = measure_SSI(AA{ig}, (e_ZW{ig}(act_factor_ind, :))');
        end
        avg_ssi = mean(ssi);

        resultfile = ['result/NGFA_BASS_sim3_n',num2str(nsamp),'_seed',num2str(seed), '_',num2str(rep),'.mat'];
        mk_prediction;
        % save(resultfile, 'e_Z', 'W', 'V', 'e_ZW', 'prec_e', 'prec_v', 'prec_w', 'act_factor_ind', 'ord', 'ssi', 'avg_ssi','activeK', 'training_error','test_error');
        disp(['nsamp ',num2str(nsamp),' / seed ', num2str(seed),' / run ', num2str(rep), ' / avg_ssi: ',num2str(avg_ssi),' finished.']);
        %%
        clear iter;
    end
    
end
end
