function X=lognormal(N,mu,sigma)

x=randn(N,1);

m=log(mu^2/sqrt(sigma^2+mu^2)); % mean of the associated normal dist.

s=sqrt(log((sigma/mu)^2+1)); % std dev of the associated normal dist.

X=exp(m+x*s);
