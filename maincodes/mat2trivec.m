% extract the upper/lower trigangular matrix of a square matrix andvers
function vectri_A = mat2trivec(A)
% without diagonal elements

n=size(A,1);
% M=n*(n-1)/2;
% N=n^2;
onetri=  tril(ones(n))';
tril_index = (onetri(:)~=1);
vec_A = A(:);
vectri_A = vec_A(tril_index)';

