% BIC for glasso
function bic_val= BIC_glasso( theta, S, num_sample, option)

p= size(theta,1);
n= num_sample;
if isempty(option)
    cmplx= log(n)*p^2/ n;
else
    switch option
        case{'r'}
            cmplx= log(n)*p^2/ n;
        case{'p'}
             % ref: StARS by Han liu, etc.(2010)
            m= p^2-sum(reshape(tril(theta), p^2, 1)==0);
             cmplx= log(n)*(2*m-p)/n;
    end
end
fitting =log(det(theta))-trace(S*theta);% -p*log(2*pi)
bic_val = -fitting+cmplx;
