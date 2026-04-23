clear
close all

%%%%%%%%%% CODE FOR ORR-SOMMERFELD EXAMPLE (FIGURE 10.8) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set parameters

R = 5772.22;            % Reynolds number
omega = 0.264002;       % perturbation frequency
N = 200;                % size of truncation
nquad = 50:50:650;      % vector for quadrature nodes for convergence test

l = 42;                 % number of left probes
rr = 42;                % number of right probes
p = 1;                  % number of moments used
mreal = 37;             % number of eigenvalues in the contour
wght = @(bb) (bb<201);  % 1./(abs(bb)).^0; % decay for the Gaussian process

rng(0);
F = diag(wght(1:N))*randn(N,l);   % sketching matrices
G = diag(wght(1:N))*randn(N,rr);  % sketching matrices

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

%% Contour integrals

cntr = 0.7+0.4i; % centre of contour
LL = 0.65; % radius of contour
Z0 = cntr;

contour = @(t) cntr + LL*cos(t)+LL*1i*sin(t);
jacobian = @(t) -LL*sin(t)+LL*1i*cos(t);

ct = 1;
Bmat = cell(2,length(nquad)); Cmat = cell(1,length(nquad));
Er = zeros(1,length(nquad));

for n = nquad % number of quadrature points
    tpts = linspace(0,2*pi,n+1)+2*pi/(2*n); tpts(end)=[];  % trap rule
    poles = contour(tpts);
    quadwts = (2*pi)/n * ones(1,n);
    residues = quadwts(:).*jacobian(tpts(:));
    [Bmat{1,ct},Bmat{2,ct},Cmat{ct}] = hankel_cont(p,l,rr,F,G,poles,residues,Z0,D2,I,Q,R,N,omega);

    % compute eigenvalues, eigenvectors, and residuals
    [V0, S0, W0] = svd(Bmat{1,ct}, 0);
    V0 = V0(:,1:mreal);  S0 = S0(1:mreal,1:mreal); W0 = W0(:,1:mreal);
    [VV, lam] = eig(V0'*Bmat{2,ct}*W0,S0,'vector');
    V=Cmat{ct}*W0*VV;
    
    lam = lam + Z0;

    Er(ct) = max(nlevp_residual(S1,D,D2,I,Q,R,N,omega,lam,V));
    ct = ct+1;

end


%% Plots

figure
plot(real(lam),imag(lam),'.k','markersize',16)
hold on
plot(real(poles),imag(poles),':k','linewidth',2)
ax=gca; ax.FontSize=18;
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',24)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',24)
legend({'eigenvalues','contour'},'fontsize',24,'interpreter','latex','location','northeast')

box on
set(gca,'layer','top');
hold on
axis equal


[V0, S0, W0] = svd(Bmat{1,end}, 0);
figure
semilogy(diag(S0)./S0(1,1),'k.--','linewidth',2,'markersize',30)
grid on
hold on
plot([mreal,mreal],[10^(-15),1],':k','linewidth',2)
ax=gca; ax.FontSize=18;
xlim([0,rr+1])
xlabel('$j$','interpreter','latex','fontsize',24)
title('$\sigma_j(\tilde{B}_0)/\sigma_1(\tilde{B}_0)$','interpreter','latex','fontsize',24)


figure
semilogy(nquad,Er,'k.-','linewidth',2,'markersize',30)
grid on
ax=gca; ax.FontSize=18;
xlabel('$n$ (number of quadrature points)','interpreter','latex','fontsize',24)
title('Max Eigenpair Residual','interpreter','latex','fontsize',24)
xlim([0,700])


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

function [B0,B1,C] = hankel_cont(p,l,rr,F,G,poles,residues,Z0,D2,I,Q,R,N,omega)

    A = cell({2*p});
    for jj = 0:2*p-1
        A{jj+1} = zeros(N-100,rr);
    end
    G=Q(1:(N-100),1:(N-100))*G;
    n = length(poles(:));
    for j = 1:n
    
        z = poles(j);
        B = -D2 + z^2*I;
        A0 = B*B/R + 1i*leg_multmat(N,chebfun(@(x) z*(1-x^2)-omega),0.5)*B  - 2i*z*I;
        B0 = [ones(1,N);  (-1).^(1:N);  B(1:end-2,:)];   
        L = A0*Q;
        [ii,~,~] = find(L(:,1));
        L = L(:,1:end-(max(ii)));
        L = [zeros(2,size(L,2));L(1:end-2,:)];
        L = B0\L;
        L = full(L(1:(N-100),1:(N-100)));
        % L=full(L);   
    
        % L=Q(:,1:(N-100))\L;

       
        U = (L\G);
    
       
        for jj = 0:2*p-1
            A{jj+1} = A{jj+1} + residues(j)*(poles(j)-Z0)^jj*U/2i/pi;
        end
    end
    
    B0 = zeros(p*l,p*rr); B1 = B0;
    C = [];
    for aa = 1:p
        if aa == 1
            for bb = 1:p
                C = [C,A{bb}];
            end
        end

        for bb = 1:p
            B0((l*(aa-1)+1):aa*l,(rr*(bb-1)+1):bb*rr) = F'*A{aa+bb-1};
            B1((l*(aa-1)+1):aa*l,(rr*(bb-1)+1):bb*rr) = F'*A{aa+bb};
        end
    end
end

function [RES] = nlevp_residual(S1,D,D2,I,Q,R,N,omega,lam,V)

RES = 0*lam;

for jj=1:length(lam)
    z = lam(jj);
    B = -D2 + z^2*I;
    A0 = B*B/R + 1i*leg_multmat(N,chebfun(@(x) z*(1-x^2)-omega),0.5)*B  - 2i*z*I;
    B3 = [ones(1,N);  (-1).^(1:N);  B(1:end-2,:)];   
    L = A0*Q;
    [ii,~,~] = find(L(:,1));
    L = L(:,1:end-(max(ii)));
    L = [zeros(2,size(L,2));L(1:end-2,:)];


    L = B3\L;
    L = L(:,1:(N-100));
    
    L=full(L);
    L = [S1;S1*D]*L;
    RES(jj) = norm(L*V(:,jj))/norm(V(:,jj));

end





end

