% estimate the model selection criteria
function mscvalue = cal_msc(hmm,crit_type) % or 'MMDL'

loglik = hmm.loglik;
invsigma = hmm.invsigma;
state_freq = hmm.state_freq;
n = sum(state_freq);
K = size(invsigma,3);
p = size(invsigma,1);
df_vec = zeros(K,1);
for i=1:K
    df_vec(i) = p^2-sum(reshape(tril(invsigma(:,:,i)), p^2, 1)==0);
end

switch crit_type
    case 'BIC'
        mscvalue = -loglik+ 1/2* log(n)*(K*(K-1) + sum(df_vec));
    case 'MMDL'
        mscvalue = -loglik + 1/2 * log(n)*K*(K-1) + 1/2* log(abs(state_freq))'*df_vec;
end