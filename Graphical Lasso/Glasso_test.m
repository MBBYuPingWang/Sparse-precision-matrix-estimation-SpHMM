% test the Graphical lasso

% model selection with AIC or BIC
clear,
 load('real_data.mat');
pop= subj_image_time{753};
S= emp_cov(pop,'row');
LambdaList =[0.01,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.4,0.5,1];
L= length(LambdaList);
[num_sample,p]= size(pop);
for i=1:L
 lambda= LambdaList(i);
[w, theta, iter, avgTol, hasError] = GraphicalLasso(pop, lambda);
%bic_val(i)= BIC_glasso(theta, S, num_sample,'p')-(log(det(theta))-trace(S*theta));
loglik(i)= log(det(theta))-trace(S*theta);
d_free(i)= 2*(sum(reshape(tril(theta), p^2, 1)~=0))/num_sample;
figure,
imagesc(w),
%colorbar
end
figure, 
subplot(1,2,1),
plot(-loglik,'r-');
hold on,
plot(d_free,'b*');
plot(-loglik+d_free);
subplot(1,2,2)
hold on,
plot(d_free/2*log(num_sample),'b*');
bic= -loglik+d_free/2*log(num_sample);
plot(-loglik,'r-');
plot(bic);
%% test the for a list of lambda
% lambdaList= 0.1:0.1:0.5;
% tic,
% [wList, thetaList, lambdaList, errors] = GraphicalLassoPath(pop, lambdaList); 
% time= toc;
% L = length(lambdaList);
% figure,
% for i=1:L
% subplot(1,L,i),
% imagesc(wList(:,:,1)),
% colorbar
% end
%  approximate, verbose, penalDiag, tolThreshold, maxIter)
