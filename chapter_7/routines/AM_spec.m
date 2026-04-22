function [sig] = AM_spec(p,q,lambda,er)
% Computes the union spectra of AM operator valid for alpha = p/q and lambda in
% [0,1]. WARNING: as coded it will only be valid for lambda at most 1
% (this is ok by Aubry duality)

if nargin < 3
    lambda=1;
    er=0;
elseif nargin <4
    er=0;
end

r = gcd(p,q);
p = round(p/r);
q = round(q/r);

if q == 2 % to take care of degenerate case
    q = 2*q;
    p = 2*p;
    theta1 = 0;
    theta3 = pi/(q/2);
else
    theta1 = 0;
    theta3 = pi/q;
end

V1 = 2*lambda*cos(2*pi*(p/q)*(1:q)+theta1);
V3 = 2*lambda*cos(2*pi*(p/q)*(1:q)+theta3);


sigma1 = sig_periodic(V1);
sigma3 = sig_periodic(V3);


sig=[sigma1;sigma3];
[~,I]=sort(sig(:,1),'ascend');
sig=sig(I,:);
sig(:,1)=sig(:,1)-er;
sig(:,2)=sig(:,2)+er;

gaps = sig(2:end,1)-sig(1:end-1,2);
tol = 100*10^(-15);
spurious = find(gaps<tol);

while ~isempty(spurious)
    sig(spurious(1),2)=max(sig(spurious(1),2),sig(spurious(1)+1,2));
    sig(spurious(1)+1,:)=[];
    gaps = sig(2:end,1)-sig(1:end-1,2);
    tol = 100*10^(-15);
    spurious = find(gaps<tol);

end