% Simulate A HMM based fMRI time courses 
%  use the HMM trained from the task related data

% load the task data
addpath('C:\Users\admin\OneDrive - Tulane University\HMM project\Refined codes');

% aligned_emoid_fmri_power264 = emoid_fmri_power264;
% for i = 1: length(emoid_fmri_power264)    
% subji_TC = emoid_fmri_power264(1).img_time_serie;
% new_subji_TC = subji_TC(Ori_index,:);
% aligned_emoid_fmri_power264(i).img_time_serie = new_subji_TC;
%  %figure, imagesc(corr(new_subj1_TC'));
% end
% save('aligned_emoid_fmri_power264.mat', 'aligned_emoid_fmri_power264', '-v7.3');
% extract part of the data
chosen_module = [{'Sensory/somatomotor Hand'};...
    {'Cingulo-opercular Task Control'};...
    {'Auditory'};...
    {'Default mode'};...
    {'Visual'}; ...
    {'Fronto-parietal Task Control'};...
    {'Salience'}];
chosen_roi = [{1:6}; {36:42}; {58:62}; {65:76}; {130:136}; {164:172}; {191:194}];
n_subj = 200;
n_T = 210;
num_module = length(chosen_module);
roi_index = [];
for module = 1:num_module
    roi_index = [roi_index, chosen_roi{module, 1}];
end

n_roi = length(roi_index);
chosen_subj_TC = zeros(n_subj, n_roi, n_T);

for i = 1: n_subj
    chosen_subj_TC(i,:,:) = aligned_emoid_fmri_power264(i).img_time_serie(roi_index,:);
end

chosen_subj_TC = normalize(chosen_subj_TC, 3);
save('chosen_subj_TC.mat', 'chosen_subj_TC');

% train a HMM on the task data
addpath(genpath('C:\Users\admin\OneDrive - Tulane University\HMM project\Refined codes'));
addpath(genpath('C:\Users\admin\Documents\MATLAB\Utilities\GroupICATv4.0b'))
% use sliding window initialize
win_size = 30;
method = 'L1';  % None
sample_size = 200;
[outFiles] = getSL(permute(chosen_subj_TC(1:sample_size,:,:), [1 3 2]), win_size, method);
% for i = 1:length(outFiles)
%     pnc_SL_data(i).SL = outFiles{1,i};
% end
dmethod = 'cityblock';
num_cluster = 3;
[Index, Centroids] = kmeans_SL(outFiles, num_cluster, dmethod);
n_roi = size(Centroids, 3);
Sigma0 = permute(Centroids, [2 3 1]);
% other initialize parameter
mu0 = randn(n_roi, num_cluster);
prior0 = ones(num_cluster,1)/num_cluster;
transmat0 = 0.5/(num_cluster-1)* ones(num_cluster);
for i=1:num_cluster
    transmat0(i,i) = 1/2;
end
mixmat0 = ones(num_cluster,1)/ num_cluster;
[LL, prior1, transmat1, mu1, Sigma1, mixmat1] = ...
    spmhmm_em(permute(chosen_subj_TC, [2 3 1]), prior0, transmat0, mu0, Sigma0, mixmat0, 'max_iter', 100);
GT_HMM_model = summary_as_hmm(LL, prior1, transmat1, mu1, Sigma1, mixmat1);
 plot_covariance_matrix(Sigma0);
 plot_covariance_matrix(Sigma1);
Initial_model = summary_as_hmm([],prior0, transmat0, mu0, Sigma0, mixmat0);
save('GT_HMM_model.mat', 'GT_HMM_model', 'Initial_model');
%% use the trained HMM to do resampling/ simulation data analyzing
addpath(genpath('C:\Users\admin\Documents\MATLAB\Utilities\HMMall'));

numex = 30;
nT = 200;
simu_obs = zeros(50, nT, numex);
simu_hidden = zeros(nT, numex);
for i = 1: numex
Sigma1 = GT_HMM_model.sigma + randn(size(GT_HMM_model.sigma))* 0.05;
[simu_obsi, simu_hiddeni] = mhmm_sample(nT, 1, ...
    GT_HMM_model.prior, GT_HMM_model.transmat, GT_HMM_model.mu, Sigma1, GT_HMM_model.mixmat);
simu_obs(:,:,i) = simu_obsi;
simu_hidden(:, i) = simu_hiddeni;
end
%%
save('simu_taskdata5.mat', 'simu_obs', 'simu_hidden'); 

%%%%%%%%% train sliding window and HMM on simu_obs seperately %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% train Sliding window
win_size = 30;
method = 'L1';
simu_obs = normalize(simu_obs,2);
[outFiles] = getSL(permute(simu_obs, [3 2 1]), win_size, method);
dmethod = 'cityblock';
num_cluster = 3;
[Index, Centroids, ] = kmeans_SL(outFiles40, num_cluster, dmethod);
plot_covariance_matrix(Centroids);
save('simutask_SLW_results', 'outFiles','Index', 'Centroids');

% train HMM
rng(3);
num_cluster = 3;
sample_size = 1;
win_size = 30;
cov_type = 'full';
[mu0, Sigma0] = mixgauss_init(num_cluster, simu_obs, cov_type);
%[mu0, Sigma0] = init_mu_sigma(simu_obs, outFiles30, num_cluster, sample_size, win_size);
n_roi = size(mu0,1);
% other initialize parameter
%
% mu0 = randn(n_roi, num_cluster);
prior0 = ones(num_cluster,1)/num_cluster;
transmat0 = 0.5/(num_cluster-1)* ones(num_cluster);
for i=1:num_cluster
    transmat0(i,i) = 1/2;
end
mixmat0 = ones(num_cluster,1)/ num_cluster;
%plot_covariance_matrix(Sigma0);
[LL, prior1, transmat1, mu1, Sigma1, mixmat1] = ...
    spmhmm_em(simu_obs, prior0, transmat0, mu0, Sigma0, mixmat0, 'max_iter', 20);

plot_covariance_matrix(Sigma1);

subj_path = decode_viterbi(simu_obs, mu1, Sigma1, prior1, transmat1);

%simu_obs= normalize(simu_obs, 2);
%   The lambda = 0.5  for the result in TMI paper 2n simulation %%%%%%%% 
trained_HMM_model = summary_as_hmm(LL, prior1, transmat1, mu1, Sigma1, mixmat1);

% switch the index
final_model = align_states(trained_HMM_model, [3 1 2]);

% decode the hidden states
subj_path = decode_viterbi(simu_obs, final_model.mu, final_model.sigma, final_model.prior, final_model.transmat);
plot_hidden_state(subj_path(1,:));
save('final_model.mat', 'final_model', 'subj_path');