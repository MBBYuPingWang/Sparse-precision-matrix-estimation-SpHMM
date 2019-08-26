% check the convergence of the EM algorithms
function converged = esti_converge(sigma,sigma2,err_bd, ...
    postmix, mixmin, previous_loglik, loglik, log_thresh)

[a, b, c]= size(sigma);
sigma_vec = reshape(sigma,[a*b*c 1]);
sigma2_vec = reshape(sigma2, [a*b*c 1]);

converged = min(postmix)/sum(postmix) < mixmin;
if converged
    fprintf('Converged: States collapsed...\n ');
else
    abs_vec = abs(sigma_vec-sigma2_vec);
    converged = max(abs_vec./(1+abs(sigma2_vec))) < err_bd;
    if converged
       fprintf('Converged: Covariances converged...\n ');
    else
       converged = em_converged(loglik, previous_loglik, log_thresh);
       if converged
           fprintf('Converged: Loglik converged...\n ');
       end
    end
end
