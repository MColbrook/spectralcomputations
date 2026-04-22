function [siglam_k,E] = fib_siglam(k, lam, symbol, d)

% function siglam_k = fib_siglam(k, lam, symbol)
%      k = level of iteration (default k = 4)
%    lam = scaling parameter (default lam = 1)
% symbol = starting symbol for substitution ('A' or 'B')
%          (default symbol = 'A')
% d is an optional input specifying number of digits for extended precision
%
% siglam_k = spectrum of level-k periodic approximation to the
%            discrete Schroedinger operator w/Fibonacci potential
% The rows of siglam_k contain the intervals that comprise the spectrum 

if nargin<1, k=4; end
if nargin<2, lam=1; end
if nargin<3, symbol='a'; end
if nargin>3
    mp.Digits(d);
    lam = mp(lam);
    tol = 100*10^(-d+1);
else
    tol = 100*10^(-15);
end
if (symbol=='a')||(symbol=='A')||(symbol==1)
    sym=0;
else
    sym=1;
end

if sym==0
    V = 1;
else
    V = 0;
end

kmax = k;
for k=2:kmax
    Vnew = zeros(fibnum(k-sym),1);
    indx = 1;
    for j=1:length(V)
        if V(j)==1, Vnew(indx:indx+1) = [1;0]; indx = indx+2; 
        else, Vnew(indx) = 1; indx = indx+1;
        end
    end
    V = Vnew;
end



N = fibnum(kmax-sym);
p = zeros(N,1);
if rem(N,2)==0
    p(1:2:end) = (1:round(N/2))';
    p(2:2:end) = (N:-1:round(N/2)+1)';
else
    p(1:2:end) = (N:-1:ceil(N/2))';
    p(2:2:end) = (1:floor(N/2))';
end

a = ones(N,1);
if nargin>3
    V = mp(V);
    a = mp(a);
end

if (kmax==1)&&(sym==0)       % A, k=1 
    siglam_k = [-2+lam 2+lam];
elseif (kmax==2)&&(sym==0)   % A, k=2
    siglam_k = [(lam-sqrt(lam^2+16))/2 0; lam (lam+sqrt(lam^2+16))/2];
elseif (kmax==1)&&(sym==1)   % B 
    siglam_k = [-2 2];
elseif (kmax==2)&&(sym==1)   % B
    siglam_k = [-2+lam 2+lam];
elseif (kmax==3)&&(sym==1)   % B
    siglam_k = [(lam-sqrt(lam^2+16))/2 0; lam (lam+sqrt(lam^2+16))/2];
else
    J1 = spdiags([a lam*V a], -1:1, N, N);
    J2 = J1;
    J1(1,N) = (-1)^N;
    J1(N,1) = (-1)^N;
    J2(1,N) = (-1)^(N+1);
    J2(N,1) = (-1)^(N+1);
    
    J1 = J1(p,p);
    J2 = J2(p,p);

    if nargin>3
        ew1 = eig(mp(J1));
        ew2 = eig(mp(J2));
    else
        ew1 = eig(J1);
        ew2 = eig(J2);
    end

    ends = sort([ew1;ew2]);
    sig = reshape(ends',2,N)';
    E = max(abs(sig(:,2)-sig(:,1)));
    gaps = sig(2:end,1)-sig(1:end-1,2);

    spurious = find(gaps<tol);
    NN = N - length(spurious);

    siglam_k = zeros(NN,2);
    if nargin>3
        siglam_k = mp(siglam_k); 
    end
    
    r1 = setdiff(1:N,spurious+1);
    r2 = setdiff(1:N,spurious);
    siglam_k(:,1) = sig(r1,1);
    siglam_k(:,2) = sig(r2,2);
end

end
