function [cent , radx , fnet] = PreImageStrSlit (Lc , Lk , ratio , n , tloerance  , Maxiter )
%
%
%
%%
t         = (0:2*pi/n:2*pi-2*pi/n).'; 
m         =  length(Lc);
cent      =  Lc;
radx      = (1-0.5.*ratio).*Lk;
rady      =  ratio.*radx;
%%

%%
err = inf;
itr = 0;
while (err>tloerance)
% for itr = 1:itr_total
itr  =itr+1;  
%%
for k=1:m
    et(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*(+radx(k).*cos(t)-i*rady(k).*sin(t));
    etp(1+(k-1)*n:k*n,1)   =          0.5.*(-radx(k).*sin(t)-i*rady(k).*cos(t));    
end
%%


%%
A       =  ones(size(et)); %% the function A in the gen. Neumann kernel
gam     =  imag(et);
%%
[mun , h ]  =  fbie(et,etp,A,gam,n,5,[],1e-13,200);
%%
fnet         =  gam+h+i.*mun;
wn           =  et-i.*fnet;
for k=1:m
    centk(k,1)  =  (max(real(wn((k-1)*n+1:k*n,1)))+min(real(wn((k-1)*n+1:k*n,1))))/2+...
                i.*(max(imag(wn((k-1)*n+1:k*n,1)))+min(imag(wn((k-1)*n+1:k*n,1))))/2; 
    radk(k,1)   =   max(real(wn((k-1)*n+1:k*n,1)))-min(real(wn((k-1)*n+1:k*n,1)));    
end
cent  =  cent-(centk-Lc);
radx  =  radx-(radk-Lk) ;
rady  =  ratio.*radx;
err   =  norm(centk-Lc,inf)+norm(radk-Lk,inf);
% [itr err]
error (itr,1) = err;
itrk  (itr,1) = itr;
%%
if itr>=Maxiter
    'No convergence after Maximunm number of iterations'
    break;
end
end
%%
end

