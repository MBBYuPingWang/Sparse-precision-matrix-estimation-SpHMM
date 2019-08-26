function mu = init_mu(data, index)

% data: nC* num_indx
% index : numindx *1
nC = size(data,1);
K = length(unique(index));

mu  =zeros(nC, K);

for i  =1:K
    mu(:,i) = mean(data(:, index==i), 2);
end
