function [loglik, exp_num_trans, exp_num_visits1, postmix, m, ip, op] = ...
    ess_hmm(prior, transmat, mixmat, mu, Sigma, data)
% ESS_MHMM Compute the Expected Sufficient Statistics for a MOG Hidden Markov Model.
%
% Outputs:
% exp_num_trans(i,j)   = sum_l sum_{t=2}^T Pr(Q(t-1) = i, Q(t) = j| Obs(l))
%  % sum_r sum_n=2^N xi(i,j)
% exp_num_visits1(i)   = sum_l Pr(Q(1)=i | Obs(l))   % sum_r gamma(z_1k)
%
% Let w(i,k,t,l) = P(Q(t)=i, M(t)=k | Obs(l))  % gamma(z^{(r)}_nk)
% where Obs(l) = Obs(:,:,l) = O_1 .. O_T for sequence l
% Then
% postmix(i,k) = sum_l sum_t w(i,k,t,l) (posterior mixing weights/
% responsibilities)
% m(:,i,k)   = sum_l sum_t w(i,k,t,l) * Obs(:,t,l)
% ip(i,k) = sum_l sum_t w(i,k,t,l) * Obs(:,t,l)' * Obs(:,t,l)
% op(:,:,i,k) = sum_l sum_t w(i,k,t,l) * Obs(:,t,l) * Obs(:,t,l)'


verbose = 1;

%[O T numex] = size(data);
numex = length(data);  % num of obs
O = size(data{1},1);   % num of time
Q = length(prior);    % num of state
M = size(mixmat,2);  % num of mix
exp_num_trans = zeros(Q,Q);
exp_num_visits1 = zeros(Q,1);
postmix = zeros(Q,M);
m = zeros(O,Q,M);
op = zeros(O,O,Q,M);
ip = zeros(Q,M);

mix = (M>1);

loglik = 0;
% if verbose, fprintf(1, 'forwards-backwards example # '); end
for ex=1:numex
  %  if verbose, fprintf(1, '%d ', ex); end
    %obs = data(:,:,ex);
    obs = squeeze(data{ex});
    T = size(obs,1);
    if mix
        [B, B2] = mixgauss_prob(obs, mu, Sigma, mixmat);
        [alpha, beta, gamma,  current_loglik, xi_summed, gamma2] = ...
            fwdback(prior, transmat, B, 'obslik2', B2, 'mixmat', mixmat);
    else
        B = mixgauss_prob(obs, mu, Sigma);
        [alpha, beta, gamma,  current_loglik, xi_summed] = fwdback(prior, transmat, B);
    end
    loglik = loglik +  current_loglik;
 %   if verbose, fprintf(1, 'll at ex %d = %f\n', ex, loglik); end
    
    exp_num_trans = exp_num_trans + xi_summed; % sum(xi,3);
    exp_num_visits1 = exp_num_visits1 + gamma(:,1);
    
    if mix
        postmix = postmix + sum(gamma2,3);
    else
        postmix = postmix + sum(gamma,2);
        gamma2 = reshape(gamma, [Q 1 T]); % gamma2(i,m,t) = gamma(i,t)
    end
    for i=1:Q
        for k=1:M
            w = reshape(gamma2(i,k,:), [1 T]); % w(t) = w(i,k,t,l)
            wobs = obs .* repmat(w, [O 1]); % wobs(:,t) = w(t) * obs(:,t)
            m(:,i,k) = m(:,i,k) + sum(wobs, 2); % m(:) = sum_t w(t) obs(:,t)
            op(:,:,i,k) = op(:,:,i,k) + wobs * obs'; % op(:,:) = sum_t w(t) * obs(:,t) * obs(:,t)'
            ip(i,k) = ip(i,k) + sum(sum(wobs .* obs, 2)); % ip = sum_t w(t) * obs(:,t)' * obs(:,t)
        end
    end
end
% if verbose, fprintf(1, '\n'); end
