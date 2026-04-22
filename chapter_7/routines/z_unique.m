 function z = z_unique(z,tol)

% function z = z_unique(z,tol)
%
% sorts z by real and then imaginary parts, and only returns entries that
% are unique up to a tolerance of O(tol).

 if nargin<3, tol=1e-12; end

 z = sort_ri(z,1e-12);

 n = length(z);
 j = 2; zval = z(1);
 indx = NaN*ones(n,1);
 indx(1) = 1; k = 2;
 while j<= n
    if abs(z(j)-zval)>tol       % found a new unique value
       zval = z(j); 
       indx(k) = j;
       k = k+1;
    end
    j=j+1;
 end 
 z = z(indx(1:k-1));
