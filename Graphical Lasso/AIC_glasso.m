% AIC for graphical Lasso
function aic_val= AIC_glasso(theta, S, num_sample, option)

p= size(theta,1);
n= num_sample;
if isempty(option)
    cmplx= 2*p^2/ n;
else
    switch option
        case{'r'}
            cmplx= 2*p^2/ n;
        case{'p'}
            % ref: StARS by Han liu, etc.(2010)
            m= p^2-sum(reshape(tril(theta), p^2, 1)==0);
            cmplx= 2*m/n;
    end
end
fitting = log(det(theta))-trace(S*theta); %-p*log(2*pi)
aic_val = -fitting+cmplx;
