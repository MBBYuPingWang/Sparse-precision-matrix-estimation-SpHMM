function A = vec2mat(vec_A, nc)
% Transform the vector to square matrix
index = (triu(ones(nc))==0);
index = index(:);
A = zeros(nc);
A(index) = vec_A;
A = reshape(A, [nc, nc]);
A = A + A' + eye(nc);
end