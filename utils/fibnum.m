 function Fn = fibnum(n)

 if n>65, fprintf('WARNING: numerical errors likely in Fibonacci number.\n'); end

 c1 = (1+1/sqrt(5))/2;
 c2 = 1-c1;

 phi = (1+sqrt(5))/2;
 alf = (1-sqrt(5))/2;

 Fn = round(c1*phi.^n + c2*alf.^n);
