function veccov = bigmat2trivec(allcov, nc)
shape = size(allcov);
if shape(1) ~= shape(2)
nW = shape(1);
allcov = reshape(allcov, [nW, nc, nc]);
veccov = zeros(nW, nc*(nc-1)/2);
for i = 1: nW
 veccov(i,:) =  mat2trivec(squeeze(allcov(i,:,:)));  
end

else
    allcov = permute(allcov, [3, 1, 2]);
    nW = size(allcov, 1);
allcov = reshape(allcov, [nW, nc, nc]);
veccov = zeros(nW, nc*(nc-1)/2);
for i = 1: nW
 veccov(i,:) =  mat2trivec(squeeze(allcov(i,:,:)));  
end
end