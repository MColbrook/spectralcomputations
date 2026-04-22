function cap = spec_capacity(sigma)
% Input: sigma, a finite collection of intervals. First column is left
% endpoints, second column is right endpoints.

% nn = ceil(log2(size(sigma,1)))+12;
nn = 8;
warning('off','all')
cap=CAP([(sigma(:,2)/2+sigma(:,1)/2)';(sigma(:,2)-sigma(:,1))'],2^nn);
warning('on','all')

end

function [tau] = CAP(IN,n)
% input
% IN: 1st row is centres, 2nd row is lengths

%%
Lc = IN(1,:)'; %the centers of the intervals
Lk = IN(2,:)'; %the length of the intervals

%%
ratio     =  0.5;
t         = (0:2*pi/n:2*pi-2*pi/n).';
m         =  length(Lc);
%%
% Lc is a vector of the midle points of the segments
% Lk is a vector of the length of the segments
[cent , radx , ~] = PreImageStrSlit (Lc , Lk , ratio , n , 1e-10  , 1000 );
rady        =   ratio.*radx;
for k=1:m
    et(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*(+radx(k).*cos(t)-1i*rady(k).*sin(t));
    etp(1+(k-1)*n:k*n,1)   =          0.5.*(-radx(k).*sin(t)-1i*rady(k).*cos(t));    
end

tau = transfinite(et, etp, cent);
end