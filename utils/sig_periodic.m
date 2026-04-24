function [sigma] = sig_periodic(V)
% computes spectrum of periodic discrete Schordinger operator with
% potential V

N = length(V);
p = zeros(N,1);
V = V(:);
if rem(N,2)==0
    p(1:2:end) = (1:round(N/2))';
    p(2:2:end) = (N:-1:round(N/2)+1)';
else
    p(1:2:end) = (N:-1:ceil(N/2))';
    p(2:2:end) = (1:floor(N/2))';
end

a = ones(N,1);

J1 = spdiags([a V a], -1:1, N, N);
J2 = J1;
J1(1,N) = (-1)^N;
J1(N,1) = (-1)^N;
J2(1,N) = (-1)^(N+1);
J2(N,1) = (-1)^(N+1);

J1 = J1(p,p);
J2 = J2(p,p);

ew1 = eig(J1);
ew2 = eig(J2);
ends = sort([ew1;ew2]);
sig = reshape(ends',2,N)';
gaps = sig(2:end,1)-sig(1:end-1,2);

tol = -100*10^(-15);

spurious = find(gaps<tol);
NN = N - length(spurious);
sigma = zeros(NN,2);
r1 = setdiff(1:N,spurious+1);
r2 = setdiff(1:N,spurious);
sigma(:,1) = sig(r1,1);
sigma(:,2) = sig(r2,2);

end