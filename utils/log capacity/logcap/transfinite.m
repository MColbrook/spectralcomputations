function [mu, hj] = transfinite(et,etp,alpha)
% mu = TRANSFINITE(et,etp,alpha) computes the transfinite diameter
% (the logarithmic capacity) of a compact set.
% 
% Notation:
%   Let eta_j : [0, 2*pi] -> Gamma_j be 2*pi-periodic and piecewise
%   twice continuously differentiable.
%   Let n be an even integer and discretize [0, 2*pi] by n nodes,
%   stored in the column vector tt.
%   (Usually: tt = 0 : 2*pi/n : 2*pi - 2*pi/n.)
% 
% Input:
%   et    = [eta_1(tt); ...; eta_L(tt)] (column vector)
%   etp   - derivative of eta (same format as et)
%   alpha = [alpha(1); ...; alpha(L)] with alpha(j) interior to
%           the j-th boundary component Gamma_j

L = length(alpha);  %% number of boundary components
n = length(et)/L;   %% number of nodes per boundary component

% Auxiliary functions gamma_j(t) = -log |eta(t)-alpha_j|
for k=1:L
    gamj(:,k) = -log(abs(et-alpha(k)));
end

% Compute the auxiliary functions h_j
A = ones(size(et)); %% the function A in the gen. Neumann kernel
for k=1:L
    [~,hjv(:,k)] = fbie(et,etp,A,gamj(:,k),n,5,[],1e-10,1000);
end

% Build and solve linear system for m_1, ..., m_L, log(mu)
for j=1:L
    for k=1:L
        hj(k,j) = sum(hjv(1+(k-1)*n:k*n,j))/n;
    end
end
% matA            = zeros(L+1,L+1);
% matA(1:L,1:L)   = hj;
matA            = hj;
matA(L+1,1:L)   = 1;
matA(1:L,L+1)   = -1;
vc_right        = zeros(L+1,1);
vc_right(L+1)   = 1;

x  = matA\vc_right;
mu = exp(x(L+1));
end