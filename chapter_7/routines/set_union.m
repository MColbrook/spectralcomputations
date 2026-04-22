function S = set_union(A,B,tol2,d)

% function S = set_union(A,B)
% 
% computes the union between two sets on the real line.
% assumes A and B are specified by intervals denoted by 
% the sorted rows of A and B.
% d is an optional input specifying number of digits for extended precision

a=1; s=1;

if nargin>3
    mp.Digits(d);
    S = mp([]);
    tol = 100*10^(-d+1);
else
    S = [];
    tol = 100*10^(-15);
end

A = [A;B];
[~,indx] = sort(A(:,1));
A = A(indx,:);
lA = size(A,1); 

while (a <= lA) 
    S(s,1) = A(a,1)-tol2;
    right = A(a,2)+tol2;
    while ((a<=lA)&&(A(a,1)-tol2<=right)) 
        right = max(right,A(a,2))+tol2;
        a=a+1; 
    end
    S(s,2) = right; s=s+1;
end
N = size(S,1);
gaps = S(2:end,1)-S(1:end-1,2);

spurious = find(gaps<tol);
if ~isempty(spurious)
    NN = N - length(spurious);
    S0 = zeros(NN,2);
    if nargin>3
        S0 = mp(S0);
    end
    r1 = setdiff(1:N,spurious+1);
    r2 = setdiff(1:N,spurious);
    S0(:,1) = S(r1,1);
    S0(:,2) = S(r2,2);
    %   fprintf('** set_union cut %d spurious intervals\n', N-NN); 
    S = S0;
end

end
