clear
close all

%%%%%%%%%% CODE FOR ORR-SOMMERFELD EXAMPLE (FIGURE 10.7) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set parameters
R = 5772.22;
omega = 0.264002;
N = 50; % truncation size

xpts = (0:0.025:1.1)-0.001;    ypts = -0.1:0.025:1; % increase resolution if wanted
zpts = kron(xpts,ones(length(ypts),1))+1i*kron(ones(1,length(xpts)),ypts(:));    zpts = zpts(:);		% complex points where we compute pseudospectra

%% Build operator pieces

N = N+100;
nn=0:N-1;
z = 1i+1;
I = eye(N);
D = leg_diffmat2(N); D2 = D*D; 
S1 = leg_normalize(N, 0.5);

bc = leg_multmat(N,chebfun(@(x) (x.^2-1).^2),0.5);
[Q,~] = qr([S1;S1*D]*bc(:,1:N-4),"econ");
Q = Q(1:N,:);
Q = diag(1./diag(S1))*Q;

%% Quartic eigenvalue problem from finite section

L10 = D2*D2/R + 1i*omega*D2;
L10 = L10*Q; coeffs{1} = L10(1:(N-100),1:(N-100));
L11 = -1i*leg_multmat(N,chebfun(@(x) (1-x^2)),0.5)*D2 -2i*I;
L11 = L11*Q; coeffs{2} = L11(1:(N-100),1:(N-100));
L12 = -2*D2/R  -1i*omega*I;
L12 = L12*Q; coeffs{3} = L12(1:(N-100),1:(N-100));
L13 = 1i*leg_multmat(N,chebfun(@(x) (1-x^2)),0.5);
L13 = L13*Q; coeffs{4} = L13(1:(N-100),1:(N-100));
L14 = I/R;
L14 = L14*Q; coeffs{5} = L14(1:(N-100),1:(N-100));
lambda = polyeig(coeffs{:});

%% Compute pseudospectra

RES=0*zpts+1;

pf = parfor_progress(length(zpts));
pfcleanup = onCleanup(@() delete(pf));

for jj=1:length(zpts)
    z = zpts(jj);
    B = -D2 + z^2*I;
    A0 = B*B/R + 1i*leg_multmat(N,chebfun(@(x) z*(1-x^2)-omega),0.5)*B  - 2i*z*I;
    B0 = [ones(1,N);  (-1).^(1:N);  B(1:end-2,:)];   
    L = A0*Q;
    [ii,~,~] = find(L(:,1));
    L = L(:,1:end-(max(ii)));
    L = [zeros(2,size(L,2));L(1:end-2,:)];


    L = B0\L;
    L = L(:,1:(N-100));
    
    L=full(L);
    L(abs(L)<10^(-18))=0;
   
    Lfs = [S1(1:(N-100),:);S1(1:(N-100),:)*D]*L(:,1:(N-100));
    L = [S1;S1*D]*L;

    a = [1:size(S1);(size(S1)+1):2*size(S1)]; a = a(:);
    L = L(a,:);
    
    II = find(vecnorm(L')<10^(-25));
    L(II,:)=[];
    RES(jj) = svds(L,1,'smallest');
    
    parfor_progress(pf);
end

 
RES=reshape(RES,length(ypts),length(xpts));

%% Plot results

v = (10.^(-16:0.25:10)); 

figure
contourf(reshape(real(zpts),length(ypts),length(xpts)),reshape(imag(zpts),length(ypts),length(xpts)),log10(real(RES)),log10(v));
cbh=colorbar;
cbh.Ticks=log10(10.^(-20:1:0));
cbh.TickLabels=["1e-20","1e-19","1e-18","1e-17","1e-16","1e-15","1e-14","1e-13","1e-12","1e-11",...
    "1e-10","1e-9","1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1"];
clim([-9,0])
colormap gray
axis tight
title('$N=50$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
box on
set(gca,'layer','top');
hold on
plot(real(lambda),imag(lambda),'.k','markersize',16) % finite section eigenvalues
set(gca,'layer','top');
axis([min(xpts(:)),max(xpts(:)),min(ypts(:)),max(ypts(:))])


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

