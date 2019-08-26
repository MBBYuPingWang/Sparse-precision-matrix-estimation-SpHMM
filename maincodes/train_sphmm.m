function [LL, prior1, transmat1, mu1, Sigma1, mixmat1]= train_sphmm(data, num_cluster)
% train the
% size of data  nC* nT * nS
% initialize with default setting
% data 
% num_cluster = 3;
parameters_of_initialize= struct( 'num_cluster', num_cluster, 'sample_size',0.3,...
    'win_size',30, 'method',[]);
time_series = permute(data,[3,2,1]);
parameters = initialize_parameters_sphmm(time_series, parameters_of_initialize);
P_init = ones(num_cluster,1)/num_cluster;
transmat = 0.5/(num_cluster-1)* ones(num_cluster);
for i=1:num_cluster
    transmat(i,i) = 1/2;
end
parameters.P_init = P_init;
parameters.transmat = transmat;

% EM iteration to fit the model
prior0 = parameters.P_init;
transmat0 = parameters.transmat;
mu0 = parameters.initMu;
Sigma0 = parameters.initcov;
mixmat0 = ones(num_cluster,1)/ num_cluster; 
[LL, prior1, transmat1, mu1, Sigma1, mixmat1] = ...
    spmhmm_em(data, prior0, transmat0, mu0, Sigma0, mixmat0, 'max_iter', 20);