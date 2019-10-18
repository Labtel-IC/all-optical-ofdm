function [y]=puls_gau(x,xd,x0,A)

y = A*exp(-1.*((x-xd)./x0).^2);