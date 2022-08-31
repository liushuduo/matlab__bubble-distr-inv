function p = lognormal_dist(R,m,s)
% LOGNORMAL_DIST value of log-normal distribution pdf at R
% variable description:
%   R:      value of R.V.
%   m:      mean of log-normal distribution
%   s:      standard variation of log-normal distribution
%   mu:     mean of log R.V.
%   sig:    standard varation of log R.V.

sig=sqrt(log(1+(s/m)^2));
mu=log(m)-0.5*sig^2;

p = exp(-(log(R)-mu).^2/(2*sig^2))./(R*sig*sqrt(2*pi));
