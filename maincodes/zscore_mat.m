function [scored_mat] = zscore_mat(ori_mat)

shape = size(ori_mat);
scored_mat = ori_mat;

for i = shape(1)
    scored_mat(i,:) = zscore(ori_mat(i,:));
end