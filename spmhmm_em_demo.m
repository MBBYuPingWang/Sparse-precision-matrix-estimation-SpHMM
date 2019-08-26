
if 0
  O = 4;
  T = 10;
  nex = 50;
  M = 2;
  Q = 3;
else
  O = 5;          %Number of coefficients in a vector 
  T = 420;         %Number of vectors in a sequence 
  nex = 3;        %Number of sequences 
  M = 1;          %Number of mixtures 
  Q = 3;          %Number of states 
end
cov_type = 'full';

data = randn(O,T,nex);

% initial guess of parameters
prior0 = normalise(rand(Q,1));

transmat0 = mk_stochastic(rand(Q,Q));

if 0
  Sigma0 = repmat(eye(O), [1 1 Q M]);
  % Initialize each mean to a random data point
  indices = randperm(T*nex);
  mu0 = reshape(data(:,indices(1:(Q*M))), [O Q M]);
  mixmat0 = mk_stochastic(rand(Q,M));
else 
  [mu0, Sigma0] = mixgauss_init(Q*M, data, cov_type);
  mu0 = reshape(mu0, [O Q M]);
  Sigma0 = reshape(Sigma0, [O O Q M]);
  mixmat0  = ones(Q,1)/Q;
  % mixmat0 = mk_stochastic(rand(Q,M));
end

plot_covariance_matrix(Sigma0);
[LL, prior1, transmat1, mu1, Sigma1, mixmat1] = ...
    spmhmm_em(data, prior0, transmat0, mu0, Sigma0, mixmat0, 'max_iter', 50);
plot_covariance_matrix(Sigma1);

loglik = mhmm_logprob(data, prior1, transmat1, mu1, Sigma1, mixmat1);


% prior0 = ones(Q,1)/ Q;
% transmat = 

%% try initialize the mu0;
% loglik = zeros(30,1);
% loglik = -Inf;
% log_vec = zeros(1, 50);
% 
%  indices = randperm(T*nex);
%   [mu0, Sigma0] = mixgauss_init(Q*M, data, cov_type);
%   mu0 = reshape(mu0, [O Q M]);
%   Sigma0 = reshape(Sigma0, [O O Q M]);
% 
%   % random ulterlize data  
% 
%   [LL, prior0, transmat0, mu0, Sigma0, mixmat0] = ...
%     spmhmm_em(data(:,:), prior0, transmat0, mu0, Sigma0, mixmat0, 'max_iter', 50);
% 
%   [LL, prior2, transmat2, mu2, Sigma2, mixmat2] = ...
%     spmhmm_em(data, prior0, transmat0, mu0, Sigma0, mixmat0, 'max_iter', 50);
% 
% figure, plot([log_vec*3/2, LL], '-*');
% plot_covariance_matrix(Sigma2, [-0.2, 0.6]);
% 
%    if loglik < LL(end)
%        Sigma = Sigma1;
%        loglik = LL(end);
%    end
