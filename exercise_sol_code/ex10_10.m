clear
close all

%% Set parameters of DDE
d1 = 2;
d2 = 1/2;
tau = 1; 
r = sqrt(23);
r1 = 6;
r2 = 3;
N = 250; % discretisation

%% Find steady state
x = chebfun('x');
cheboppref.maxIter=50;
L = chebop(@(x,u,v) [d1*tau*diff(u,2) + tau*u.*(r1-v); d2*tau*diff(v,2) + tau*v.*(-r2+u)]);
L.lbc = @(u,v)[ u-2; v-1];
L.rbc =  @(u,v)[ u; v];
L.init = [0*x; 1+x.^2];
[u,v] = L\0;

%% Build operator pieces
N = N+300; % to deal with bands, reduce by 300 later
D = leg_diffmat2(N+2); D2 = D*D;
U = leg_multmat(N+2,u,0.5);
V = leg_multmat(N+2,v,0.5);
S1 = leg_normalize(N+2, 0.5);
I = eye(N+2);

%% Take care of boundary conditions

bc = leg_multmat(N+2,chebfun(@(x) (x.^2-1)),0.5);
[Q,~] = qr(S1*bc(:,1:N),"econ");
Q = diag(1./diag(S1))*Q;

D2 = D2*Q;
U = U*Q;
V = V*Q;
I = I*Q;

c = bandwidth(D2,'lower'); c = max(c,bandwidth(U,'lower')); c = max(c,bandwidth(V,'lower')); c = max(c,bandwidth(I,'lower'));

N = N-300;
S1_fin = S1(1:(N),1:(N));
S1 = S1(1:(N+c),1:(N+c));

D2_fin = S1_fin*D2(1:N,1:N);
D2 = S1*D2(1:(N+c),1:N);
U_fin = S1_fin*U(1:N,1:N);
U = S1*U(1:(N+c),1:N);
V_fin = S1_fin*V(1:N,1:N);
V = S1*V(1:(N+c),1:N);
I_fin = S1_fin*I(1:N,1:N);
I = S1*I(1:(N+c),1:N);


%% Build pieces of pencil

LI_fin = [I_fin,0*I_fin;0*I_fin,I_fin];
LI = [I,0*I;0*I,I];

Ld_fin = tau*[d1*D2_fin,   0*I_fin;
    0*I_fin,   d2*D2_fin];
Ld = tau*[d1*D2,   0*I;
    0*I,   d2*D2];

L1_fin = tau*[r1*I_fin-V_fin,   0*I_fin;
    0*I_fin,   -r2*I_fin+U_fin];
L1 = tau*[r1*I-V,   0*I;
    0*I,   -r2*I+U];

L2a_fin = tau*[0*I_fin,   -U_fin;
    0*I_fin,   0*I_fin];
L2a = tau*[0*I,   -U;
    0*I,   0*I];

L2b_fin = tau*[0*I_fin,   0*I_fin;
    V_fin,   0*I_fin];
L2b = tau*[0*I,   0*I;
    V,   0*I];

%% Contour integral method
close all
l = 30;
rr = 30;
mreal = 24;

% contour parameters
cntr = 2;
Z0 = cntr;
LL = 3.45;
contour = @(t) cntr + LL*cos(t)+LL*1i*sin(t);
jacobian = @(t) -LL*sin(t)+LL*1i*cos(t);

% Loewner method parameters
contour2 = @(t) cntr + 1.04*(LL*cos(t)+LL*1i*sin(t));
tpts2 = linspace(0,2*pi,2*l+1); tpts2(end)=[];
tpts2 = tpts2 - pi+ 2*pi/(4*l);
tpts2 = tpts2/4.5 -pi;

zleft = contour2(tpts2(1:2:end));
zright = contour2(tpts2(2:2:end));


% random probing matrices
rng(0);
wght = @(bb) 1./(abs(bb));
F = [diag(wght(1:N))*(randn(N,l)+1i*randn(N,l))/sqrt(2);diag(wght(1:N))*(randn(N,l)+1i*randn(N,l))/sqrt(2)];
G = [diag(wght(1:N))*(randn(N,rr)+1i*randn(N,rr))/sqrt(2);diag(wght(1:N))*(randn(N,rr)+1i*randn(N,rr))/sqrt(2)];

nquad = 1000;
ct = 1;
hcn = zeros(length(nquad),4); hres = hcn;
lres = zeros(length(nquad),1); lcn = lres;
n = nquad;
tpts = linspace(0,2*pi,n+1)+2*pi/(2*n); tpts(end)=[];  % trap rule
quadpts = contour(tpts);
quadwts = (2*pi)/n * ones(1,n);
totwts = quadwts(:).*jacobian(tpts(:));

[M,Mv] = pred_prey_solve(l,rr,LI_fin,L1_fin,Ld_fin,L2a_fin,L2b_fin,r,quadpts,F,G);

hleft = cell(1,1);
hleft{1,1} = @(z) (z-Z0).^0;
hleft{2,1} = @(z) (z-Z0).^1; 
hleft{3,1} = @(z) (z-Z0).^3;

hleft{4,1} = @(z) (z-Z0).^4;
[V,zpts,hcn] = Hankel_cont(M,Mv,Z0,quadpts,totwts,hleft,hleft,mreal);

%%
[zpts,V] = pred_prey_newton(LI_fin,L1_fin,Ld_fin,L2a_fin,L2b_fin,r,zpts,V,4,V);
Res = pred_prey_residual(LI,Ld,L1,L2a,L2b,V,zpts);


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



function [M,Mv] = pred_prey_solve(l,rr,LI_fin,L1_fin,Ld_fin,L2a_fin,L2b_fin,r,quadpts,F,G)
n = length(quadpts(:));
M = zeros(l,rr,n);
Mv = zeros(size(F,1),rr,n);
for j = 1:n
    Mv(:,:,j) = (  (quadpts(j)*LI_fin - Ld_fin - L1_fin -L2a_fin*exp(-quadpts(j)*r) -L2b_fin*exp(-quadpts(j)))\G           );
    M(:,:,j) = F'*Mv(:,:,j);
end
end


function [Z0,V0] = pred_prey_newton(LI_fin,L1_fin,Ld_fin,L2a_fin,L2b_fin,r,Z0,V0,NI,W)
W = W./vecnorm(W);
V0 = V0./(dot(W,V0));
n = length(Z0(:));
warning('off','all')
for ii = 1:NI
    for j = 1:n
        v0 = (LI_fin + r*L2a_fin*exp(-Z0(j)*r) + L2b_fin*exp(-Z0(j)))*V0(:,j);
        v1 = (Z0(j)*LI_fin - Ld_fin - L1_fin -L2a_fin*exp(-Z0(j)*r) -L2b_fin*exp(-Z0(j)))\v0;
        Z0(j) = Z0(j)-1/(W(:,j)'*v1);
        V0(:,j) = v1/(W(:,j)'*v1);
    end
end
warning('on','all')
end

function RES = pred_prey_residual(LI,Ld,L1,L2a,L2b,V,zpts)
RES = 0*zpts;
r = sqrt(23);
    
for jj=1:length(zpts)
    L = zpts(jj)*LI - Ld - L1 -L2a*exp(-zpts(jj)*r) -L2b*exp(-zpts(jj));
    RES(jj) = norm(L*V(:,jj))/norm(V(:,jj)); 
end


end

function [V,lam,c] = Hankel_cont(M,Mv,Z0,quadpts,totwts,hleft,hright,mreal)

p =size(hleft,1);
n = length(quadpts);
l = size(M,1);
rr = size(M,2);

B0 = zeros(p*l,p*rr); B1 = B0;
for aa = 1:p
    for bb = 1:p
        w = zeros(1,1,n);
        w(1,1,:) = transpose(totwts(:).*hleft{aa}(quadpts(:)).*hright{bb}(quadpts(:)));
        B0((l*(aa-1)+1):aa*l,(rr*(bb-1)+1):bb*rr) = squeeze(sum(M.*w,3));

        w(1,1,:) = transpose((quadpts(:)-Z0).*totwts(:).*hleft{aa}(quadpts(:)).*hright{bb}(quadpts(:)));
        B1((l*(aa-1)+1):aa*l,(rr*(bb-1)+1):bb*rr) = squeeze(sum(M.*w,3));
    end
end

[V0, S0, W0] = svd(B0, 0);
[VV, lam] = eig(V0(:,1:mreal)'*B1*W0(:,1:mreal),S0(1:mreal,1:mreal),'vector');
[~,I] = sort(real(lam),'ascend'); VV = VV(:,I); lam = lam(I) + Z0;
c = S0(1,1)./S0(mreal,mreal);

C = zeros(size(Mv,1),p*rr);
for bb = 1:p
    w = zeros(1,1,n);
    w(1,1,:) = transpose(totwts(:).*hleft{1}(quadpts(:)).*hright{bb}(quadpts(:)));
    C(:,(rr*(bb-1)+1):bb*rr) = squeeze(sum(Mv.*w,3));
end

V=C*W0(:,1:mreal)*VV; V = V./vecnorm(V);

end

function [V,lam,c] = Loewner_cont(M,Mv,Z0,quadpts,totwts,zleft,zright,mreal)

p = length(zleft);
n = length(quadpts);

weight0 = zeros(p,p,n);
weight1 = zeros(p,p,n);
for aa = 1:p
    for bb = 1:p
        weight0(aa,bb,:) = transpose(totwts(:).*1./(quadpts(:)-zleft(aa)).*1./(quadpts(:)-zright(bb)));
        weight1(aa,bb,:) = transpose(totwts(:).*1./(quadpts(:)-zleft(aa)).*1./(quadpts(:)-zright(bb)).*(quadpts(:)-Z0));
    end
end

B0 = squeeze(sum(M.*weight0,3));
B1 = squeeze(sum(M.*weight1,3));
C = squeeze(sum(Mv.*weight0(1,:,:),3));

[V0, S0, W0] = svd(B0, 0);
[VV, lam] = eig(V0(:,1:mreal)'*B1*W0(:,1:mreal),S0(1:mreal,1:mreal),'vector');
[~,I] = sort(real(lam),'ascend'); VV = VV(:,I); lam = lam(I) + Z0;
c = S0(1,1)./S0(mreal,mreal);

V=C*W0(:,1:mreal)*VV; V = V./vecnorm(V);

end







