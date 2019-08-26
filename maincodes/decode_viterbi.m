% viterbi step for decoding the HMM
function subj_path = decode_viterbi(data, mu, Sigma, prior, transmat)
% data size : nC * nT * nS

%num_state = length(prior);
num_time = size(data,2);
num_subj = size(data,3);

subj_path= zeros(num_subj, num_time);
    for i=1: num_subj
        obslik = mixgauss_prob(data(:,:,i), mu, Sigma);
        subj_path(i, :) = viterbi_path(prior, transmat, obslik);
    end
% save(path_str,'subj_path');
