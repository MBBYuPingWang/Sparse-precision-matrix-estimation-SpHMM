% calculate emperical coviarance
function S= emp_cov(pop,option)

% option
% - vector of instance
[~,n]= size(pop);
if isempty(option)
    option='col';
else
    switch option
        case{'row'}
         mu = sum(pop)/n;
         temp = bsxfun(@minus,pop,mu);
         S= temp'*temp/n;
        case{'col'}
           mu= sum(pop,2)/n;
           temp= bsxfun(@minus,pop,mu);
           S= temp*temp'/n;
    end
end
