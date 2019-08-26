function [val] = l1norm(trade_off, tar_matrix)
L =length(trade_off);
shape = size(tar_matrix);

if length(shape)< 3 
    % l1 norm of sequence of vectors
            val = 0;
        for i=1: L
            itar_vec = tar_matrix(:,i);
            val = val + trade_off(i)*sum(abs(itar_vec(:)));
        end
elseif length(shape) ==3    
    if shape(3)== L
        val = 0;
        for i=1: L
            itar_matrix = tar_matrix(:,:,i);
            val = val + trade_off(i)*sum(abs(itar_matrix(:)));
        end
    else
        error('The dimension of input are not compatible! please check.');
    end
else
    error('The dimension of tar_matrix is at least 2. please check');
end
