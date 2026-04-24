clear
close all

%% First code up the discretisations using Legendre polynomials and plot eigenvalues

N = 1000;

% matrices
N = N + 2; % to account for the multiplication L*B later
nn=0:N-1;
D = (leg_diffmat2(N));
L = -D*leg_multmat(N,chebfun(@(x) (1+x^2)),0.5)*D;
B = -[sparse(2,N);speye(N-2,N)]+[speye(N-2,N);sparse(2,N)];
A = L*B;
N = N - 2;
A = A(:,1:N); % keep extra rows for residual later
B = B(:,1:N);

% compute eigenvalues
[V,E]=eig(full(A(1:N,:)),full(B(1:N,:)),'vector');
[E,I] = sort(E,'ascend');
V = V(:,I);

% plot results
figure
semilogy(1:N,E,'.k','markersize',12)
% hold on
% plot([N,N]*2/pi,[1,max(E)],'k--','linewidth',1)
xlabel('Eval Number','interpreter','latex','fontsize',18)
title(sprintf('$N=%d$',N),'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
grid on
ylim([1,max(E)])

exportgraphics(gcf,'leg_evals3.pdf','ContentType','vector','BackgroundColor','none','Colorspace','gray')

%% Now compute error bounds using orthogonality relation of Legendre polynomials

S1 = leg_normalize(N+2, 0.5);
Res = zeros(N,1);
for j=1:N-2
    Res(j) = norm(S1*((A-E(j)*B)*V(:,j)))/norm(S1*B*V(:,j));
end
    
% plot results
figure
semilogy(1:N,Res./E,'.k','markersize',12)
% hold on
% plot([N,N]*2/pi,[10^(-15),1],'k--','linewidth',1)
xlabel('Eval Number','interpreter','latex','fontsize',18)
title(sprintf('$N=%d$',N),'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
grid on

exportgraphics(gcf,'leg_eval3_errors3.pdf','ContentType','vector','BackgroundColor','none','Colorspace','gray')


%% Code for spectral discretisation using Legendre polynomials

function S = leg_normalize(n, lambda)
% normalisation for Legendre polynomials (used to compute L^2 residuals)
C = ones(1,n);
for ii = 1:round(2*lambda-1)
    C = C.*((1:n)+ii-1);
end
C = sqrt(((2^(1-2*lambda)*pi/(gamma(lambda)^2))./((0:n-1)+lambda)).*C);
S = spdiags(C',0,n,n);
end

function D = leg_diffmat2(n)
% Legendre differentiation matrix
D = sparse(n,n);
for jj=2:n
    nn = (jj-2):-2:0;
    D((jj-1:-2:1),jj) = 2*nn+1;
end
end

function M = leg_multmat(n, f, lambda)
% Matrix for multiplication by a chebfun f in Legendre space
a = legcoeffs(f);

% Multiplying by a scalar is easy.
if ( numel(a) == 1 )
    M = a*speye(n);
    return
end

% Prolong or truncate coefficients
if ( numel(a) < n )
    a = [a ; zeros(n - numel(a), 1)];   % Prolong
else
    a = a(1:n);                         % Truncate.
end

% Convert to C^{lam}
a = ultraS.convertmat(n, 0.5, lambda - 1) * a;

m = 2*n; 
M0 = speye(m);

d1 = [1 (2*lambda : 2*lambda + m - 2)]./ ...
    [1 (2*((lambda+1) : lambda + m - 1))];
d2 = (1:m)./(2*(lambda:lambda + m - 1));
B = [d2' zeros(m, 1) d1'];
Mx = spdiags(B,[-1 0 1], m, m);
M1 = 2*lambda*Mx;

% Construct the multiplication operator by a three-term recurrence: 
M = a(1)*M0;
M = M + a(2)*M1;
for nn = 1:length(a) - 2
    M2 = 2*(nn + lambda)/(nn + 1)*Mx*M1 - (nn + 2*lambda - 1)/(nn + 1)*M0;
    M = M + a(nn + 2)*M2;
    M0 = M1;
    M1 = M2;
    if ( abs(a(nn + 3:end)) < eps ), break, end
end
M = M(1:n, 1:n); 

end

