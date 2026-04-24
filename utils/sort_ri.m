 function [z,indx] = sort_ri(z, tol)

% function [z,indx] = sort_ri(z, tol)
%
% Sorts the list of complex numbers z by increasing real part,
% and secondarily by increasing imaginary part.
% The "precision" of the first sort key is "tol", so that
% all entries in z that have real parts that differ by tol (ish)
% are considered equal.

 if nargin<2, tol=1e-10; end

 [zr, indx1] = sort(real(z));
 z = z(indx1);
 N = length(z);

 jl = 1;                        % initial left marker
 jr = 2;                        % initial right marker
 zr = [z(1)];
 indx2 = [1:N]';
 while jr<=N
    if ~((real(z(jr))-real(z(jr-1))) <= tol)   % break off equal real values
       jj  = jl:jr-1;
       [~,indxi] = sort(imag(z(jj)));
       indx2(jj) = jl+indxi-1;
       jl=jr; 
    elseif jr==N
       jj  = jl:jr;
       [~,indxi] = sort(imag(z(jj)));
       indx2(jj) = jl+indxi-1;
    end 
    jr=jr+1;                 
 end
 z = z(indx2);
 indx = indx1(indx2);
