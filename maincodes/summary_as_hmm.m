function model = summary_as_hmm(LL, prior1, transmat1, mu1, Sigma1, mixmat1)

model  = struct('LL',LL, 'prior', prior1, 'transmat',transmat1, 'mu',mu1,...
        'sigma', Sigma1, 'mixmat',mixmat1);
    