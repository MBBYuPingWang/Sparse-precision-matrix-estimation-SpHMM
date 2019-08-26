%the objective function of glasso
function obj_val = cal_glasso(invsigma, S, lambda_list)

L = length(lambda_list);
obj_val = 0;

for l = 1:L
    invsigma_l = squeeze(invsigma(:,:,l));
    S_l = squeeze(S(:,:,l));
    obj_val = obj_val -logdet(invsigma_l) +  trace(S_l*invsigma)...
        + lambda_list(l)*sum(abs(invsigma_l(:)));
end
