function [mu0, sigma0] = init_mu_sigma(simu_obs, outFiles, num_cluster, sample_size, win_size)

n_roi = size(simu_obs, 1);
nS = ceil(length(outFiles) * sample_size);
temp = randperm(nS);
subj_indx = temp(1:nS);

sample_outFiles = outFiles(subj_indx);
sample_obs = simu_obs(:,:,subj_indx);
dmethod = 'cityblock';
[Index, SLC] = kmeans_SL(sample_outFiles, num_cluster, dmethod);
sigma0 = permute(SLC, [2 3 1]);

mu0 = init_mu(reshape(sample_obs(:,win_size+1:end,:), [n_roi, length(Index)]), Index);

