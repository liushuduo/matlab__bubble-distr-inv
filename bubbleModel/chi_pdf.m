function p=chi_pdf(x,theta,a)

p=2*x.^(a-1).*exp(-x.^2/theta)/(theta^(a/2)*gamma(a/2));