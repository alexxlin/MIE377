function [sr_exante, sr_expost] = Sharp_Ratio(x, mu, Q, mu_p, sigma_p)
    sr_exante = mu'*x/sqrt(x'*Q*x);
    sr_expost = mu_p/sigma_p;
end