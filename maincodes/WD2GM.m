function D = WD2GM(X, Y)

% The Wasserstain-2 Distance between two Gaussian distribution
A = size(X.mu,2);
B = size(Y.mu,2);
D=zeros(A,B);
for i=1:A
    cov1=X.sigma(:,:,i);
    sqrt_cov1=sqrtm(cov1);    
    for j=1:B
        cov2=Y.sigma(:,:,j);
        D(i,j)= norm(X.mu(:,i)-Y.mu(:,j))^2+...
            trace(cov1+cov2-2*sqrtm(sqrt_cov1*cov2*sqrt_cov1));
        
        D(i,j)=real(sqrt(D(i,j)));
        %D(i,j)=wasserstein_gaussian(X.mu(:,i),X.covariance(:,:,i),Y.mu(:,j),Y.covariance(:,:,j));
%         disp(D(i,j))
    end
end
D=1.0*D;
end