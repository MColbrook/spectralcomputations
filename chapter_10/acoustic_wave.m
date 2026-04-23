clear
close all

%%%%%%%%%%%%%%% CODE FOR THE ACOUSTIC WAVE EXAMPLE %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Truncation and FEM

nvec = [10,100,500];
F = zeros(4,length(nvec));
ct = 1;
lamFD = cell(length(nvec),1);
Efun = lamFD; Er = lamFD;

for n = nvec
    cfs = acoustic_wave_1d(n,1);
    [X,E] = polyeig(cfs{:});
    [~,I] = sort(abs(real(E)));
    X = flipud(X(:,I)); E = E(I);
    xx = 1/(2*n):1/n:1;
    err = zeros(1,length(E));
    for jj = 1:length(E)
        yy2 = (exp(1i*2*pi*E(jj)*xx));
        err(jj)=subspacea(X(:,jj),yy2(:)); % subspace angle of FEM vs. Plane Waves
    end
    Efun{ct} = X;
    lamFD{ct} = E;
    Er{ct} = err;

    ct = ct+1;
end


%% Figure 10.5

figure
plot(real(lamFD{1}),imag(lamFD{1}),'k+','markersize',6,'linewidth',2)
hold on
plot(real(lamFD{2}),imag(lamFD{2}),'ko','markersize',8)
plot(real(lamFD{3}),imag(lamFD{3}),'k.','markersize',15)
ax = gca; ax.FontSize = 18;
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
legend({'$n=10$','$n=100$','$n=500$'},'fontsize',16,'interpreter','latex','location','northeast')
grid on


%% Figure 10.6

figure
loglog(abs(real(lamFD{1})),Er{1},'k+','markersize',6,'linewidth',2)
hold on
loglog(abs(real(lamFD{2})),Er{2},'ko','markersize',8)
loglog(abs(real(lamFD{3})),Er{3},'k.','markersize',15)
ax = gca; ax.FontSize = 18;
xlabel('$|\mathrm{Re}(\lambda)|$','interpreter','latex','fontsize',18)
legend({'$n=10$','$n=100$','$n=500$'},'fontsize',16,'interpreter','latex','location','southeast')
title('Subspace Angle: FEM vs. Plane Waves','interpreter','latex','fontsize',18)


figure

subplot(2,2,1)
JJ = 1; n=10;
yy1=Efun{1}(:,JJ);
xx = 1/(2*n):1/n:1;
plot(xx,yy1,'k.','markersize',18)
hold on
yy2 = (exp(1i*2*pi*lamFD{1}(JJ)*xx)); yy2 = yy2(:);
c=yy2\yy1;
xx =0:0.01:1;
yy2 = (exp(1i*2*pi*lamFD{1}(JJ)*xx)); yy2 = yy2(:);
yy2 = c*yy2; % eigenfunctions only defined up to a constant so make them match as close as possible
plot(xx,yy2,'k:','linewidth',1.2)
ax = gca; ax.FontSize = 12;
xlim([0,1])
xlabel('$x/L$','interpreter','latex','fontsize',18)

subplot(2,2,2)
JJ = 3; n=10;
yy1=Efun{1}(:,JJ);
xx = 1/(2*n):1/n:1;
plot(xx,yy1,'k.','markersize',18)
hold on
yy2 = (exp(1i*2*pi*lamFD{1}(JJ)*xx)); yy2 = yy2(:);
c=yy2\yy1;
xx =0:0.01:1;
yy2 = (exp(1i*2*pi*lamFD{1}(JJ)*xx)); yy2 = yy2(:);
yy2 = c*yy2;
plot(xx,yy2,'k:','linewidth',1.2)
ax = gca; ax.FontSize = 12;
xlim([0,1])
xlabel('$x/L$','interpreter','latex','fontsize',18)

subplot(2,2,3)
JJ = 5; n=10;
yy1=Efun{1}(:,JJ);
xx = 1/(2*n):1/n:1;
plot(xx,yy1,'k.','markersize',18)
hold on
yy2 = (exp(1i*2*pi*lamFD{1}(JJ)*xx)); yy2 = yy2(:);
c=yy2\yy1;
xx =0:0.01:1;
yy2 = (exp(1i*2*pi*lamFD{1}(JJ)*xx)); yy2 = yy2(:);
yy2 = c*yy2;
plot(xx,yy2,'k:','linewidth',1.2)
ax = gca; ax.FontSize = 12;
xlim([0,1])
xlabel('$x/L$','interpreter','latex','fontsize',18)

subplot(2,2,4)
JJ = 7; n=10;
yy1=Efun{1}(:,JJ);
xx = 1/(2*n):1/n:1;
plot(xx,yy1,'k.','markersize',18)
hold on
yy2 = (exp(1i*2*pi*lamFD{1}(JJ)*xx)); yy2 = yy2(:);
c=yy2\yy1;
xx =0:0.01:1;
yy2 = (exp(1i*2*pi*lamFD{1}(JJ)*xx)); yy2 = yy2(:);
yy2 = c*yy2;
plot(xx,yy2,'k:','linewidth',1.2)
ax = gca; ax.FontSize = 12;
xlim([0,1])
xlabel('$x/L$','interpreter','latex','fontsize',18)

sgtitle('FEM Eigenfunctions','interpreter','latex','fontsize',18)


%% Compute pseudospectra of operator on half line

N = 100;
xpts=-1:0.02:1;    ypts=-1:0.02:0.25;
zpts=kron(xpts,ones(length(ypts),1))+1i*kron(ones(1,length(xpts)),ypts(:));    zpts=zpts(:);		% complex points where we compute pseudospectra

RES=0*zpts+1;
pf = parfor_progress(length(zpts));
pfcleanup = onCleanup(@() delete(pf));

for jj=1:length(zpts)
    L = mat_setup(zpts(jj),N,[]);
    RES(jj) = min(svd(L(:,1:end)));
    parfor_progress(pf);  
end

RES=reshape(RES,length(ypts),length(xpts));


%% Plot the results


v=(10.^(-20:0.2:0));
figure
contourf(reshape(real(zpts),length(ypts),length(xpts)),reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES),min(v))),log10(v),'LineColor',[1,1,1]*0,...
    'linewidth',1,'linestyle','-','ShowText','off');
cbh=colorbar;
cbh.Ticks=log10(10.^(-20:1:0));
cbh.TickLabels=["1e-20","1e-19","1e-18","1e-17","1e-16","1e-15","1e-14","1e-13","1e-12","1e-11",...
    "1e-10","1e-9","1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1"];
clim([-4,0])
colormap gray
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
title('$N=100$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
box on
set(gca,'layer','top');
% 
exportgraphics(gcf,'CHAP9_wave_N100.pdf','ContentType','vector','BackgroundColor','none','Colorspace','gray')






function [L,Z] = mat_setup(lam,N,X)

%% Set up matrices using Laguerre functions

NN = N+1;
D = zeros(NN,NN);
for jj = 1:NN
    D(1:jj-1,jj)= -1;
    D(jj,jj) = -1/2;
end

L = D*D + lam^2*eye(NN);

% Basis recombination

B = [zeros(1,NN-1);eye(NN-1,NN-1)];
for jj = 1:NN-1
     B(1,jj) = (-(jj+1/2)-1i*lam)/(0.5+1i*lam);
end


[Q,~] = qr(B,"econ");

L = L*Q;
NN = NN - 1;
L = L(:,1:NN);


if ~isempty(X)
    X=X(:);
    LL = zeros(length(X),N+1);
    LL(:,1) = exp(-X/2);
    LL(:,2) = (1 - X).*exp(-X/2);
    for jj=2:N
        k = jj - 1;
        LL(:,jj+1) = (2*k+1 - X)./(k+1).*LL(:,jj)-k/(k+1)*LL(:,jj-1);
    end
    Z = LL*Q;
else
    Z = [];
end

end


















function [coeffs] = acoustic_wave_1d(n,z)
%ACOUSTIC_WAVE_1D   Acoustic wave problem in 1 dimension.
%  [COEFFS,FUN,F] = nlevp('acoustic_wave_1d',N,Z) constructs an N-by-N
%  quadratic matrix polynomial lambda^2*M + lambda*D + K that arises
%  from the discretization of a 1D acoustic wave equation.
%  The damping matrix has the form 2*pi*i*Z^(-1)*C, where
%  C = e_n*e_n', where e_n = [0 ... 0 1]', and the scalar parameter Z is
%  the impedance (possibly complex).
%  The default values are N = 10 and Z = 1.
%  The eigenvalues lie in the upper half of the complex plane.
%  The matrices are returned in a cell array: COEFFS = {K, D, M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle K + lambda*D + lambda^2*M.
%  XCOEFFS returns the cell {1 en 1;K en' M} to exploit the low rank of D.
%  This problem has the properties pep, qep, symmetric, *-even
%  parameter-dependent, scalable, sparse, tridiagonal, banded, low-rank.

%  Reference:
%  F. Chaitin-Chatelin and M. B. van Gijzen, Analysis of parameterized
%  quadratic eigenvalue problems in computational acoustics with homotopic
%  deviation theory, Numer. Linear Algebra Appl. 13 (2006), pp. 487-512.

if nargin < 2 || isempty(z)
    z = 1; 
end
if nargin < 1 || isempty(n) 
    n = 10; 
end

h = 1/n;

e = ones(n,1);
K = spdiags([-e,2*e,-e],-1:1,n,n);
K(n,n) = 1;
K = n*K;

D = sparse(n,n,1/z,n,n);
M = speye(n); M(n,n) = 0.5; M = h*M;

coeffs = {K, 2*pi*1i*D, -(2*pi)^2*M};

end


    


