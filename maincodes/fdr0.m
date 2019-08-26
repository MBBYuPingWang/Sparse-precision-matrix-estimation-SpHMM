function y = fdr0(p, q)
% y = fdr0(p, q)
%
% to calculate whether a pvalue survive FDR corrected q
%
% p: an array of p values. (e.g. p values for each channel)
% q: desired FDR threshold (typically 0.05)
% y: an array of the same size with p with only two possible values. 0
% means this position (channel) does not survive the threshold, 1 mean it
% survives
%
% Ref:
% Genovese et al. (2002). Thresholding statistical maps in functional
% neuroimaging using the false discovery rate. Neuroimage, 15:722-786.
%
% Example:
%   y = fdr0(rand(10,1),0.5);
%
% Xu Cui
% 2016/3/14
%

pvalue = p;
y = 0 * p;

[sortedpvalue, sortedposition] = sort(pvalue);
v = length(sortedposition);
for ii=1:v
    if q*ii/v >= sortedpvalue(ii)
        y(sortedposition(ii)) = 1;
    end
end

return;