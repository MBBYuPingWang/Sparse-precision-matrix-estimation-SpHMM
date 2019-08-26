function new_model = align_states(old_model, index)

new_model = old_model;
 
new_model.mu = old_model.mu(:, index);
new_model.sigma = old_model.sigma(:,:,index);
new_model.prior =  old_model.prior(index);
new_model.transmat = old_model.transmat(index, index);
