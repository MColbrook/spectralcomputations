clear
close all

% Add SpecSolve folder from chapter 4 (constructs the rational kernels)
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'chapter_4', 'SpecSolve');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%% CODE FOR Example in Section 5.5.2 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set parameters and build matrix

n = 10^4;
beta = 1;
M = 1;
kappa = 1;

d = ones(n,1);
r = floor(n/(M+1));
d((M+1):(M+1):end) = (((1:r)+1)./(1:r)).^beta*kappa;
d = [1;d(:)];

A = spdiags(d,1,n,n); A = A+A';
A = A(1:n,1:n-1);
b = zeros(n,1);b(1)=1;

%% Compute spectral measure - use algorithms from chapter 4 (specsolve)

X=[-3:0.01:-2.6,-2.6:0.001:-2,-2:0.01:-0.04,-0.04:0.001:0];
X = unique([X,-fliplr(X)]); X= sort(X);

mu1 = infmatMeas(A,b,X,0.1,'order',2);
mu2 = infmatMeas(A,b,X,0.01,'order',2);
mu3 = infmatMeas(A,b,X,0.001,'order',2);

K = 0;
[poles,res] = rational_kernel(2,"equi");
for jj=1:2
    K = K + res(jj)./(0-poles(jj));
end
K = imag(K)/pi;

%% Plot results and visualise BIC

figure
semilogy(X,mu1,'g','LineWidth',1)
hold on
plot(X,mu2,'r','LineWidth',1.5)
plot(X,mu3,'k','LineWidth',2)
legend({'$\epsilon=0.1$','$\epsilon=0.01$','$\epsilon=0.001$'},'fontsize',16,'interpreter','latex','location','south')
xlim([min(X),max(X)])
xlabel('$x$','interpreter','latex','fontsize',18)
title('$[K_\epsilon\ast\mu_v](x)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

figure
plot(X,mu1/K*0.1,'g','LineWidth',1)
hold on
plot(X,mu2/K*0.01,'r','LineWidth',1.5)
plot(X,mu3/K*0.001,'k','LineWidth',2)
legend({'$\epsilon=0.1$','$\epsilon=0.01$','$\epsilon=0.001$'},'fontsize',16,'interpreter','latex','location','southeast')
xlim([-0.1,0.1])
xlabel('$x$','interpreter','latex','fontsize',18)
title('$[K_\epsilon\ast\mu_v](x)/K_\epsilon(0)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;


%% Functional calculus computation of pure point part

B = zeros(3,1500);
n = 2000;
H = full(A(1:n,1:n));
[V,E] = eig(H,'vector'); E = E(:); % direct diagonalisation
I = find(abs(E)<1); I = I(:);
P = V'*b(1:n); P = P(I);

pf = parfor_progress(size(B,2));
pfcleanup = onCleanup(@() delete(pf));
for T = 1:size(B,2)
    m = ceil(3*T);
    [ccn,ccw] = chebpts(m,[0,T]);
    vv = (transpose(exp(-1i*ccn*E(I)')).*P);

    n2 = 5;
    B(1,T) = vecnorm(V(1:n2,I)*vv).^2*ccw'/T;
    n2 = 10;
    B(2,T) = vecnorm(V(1:n2,I)*vv).^2*ccw'/T;
    n2 = 100;
    B(3,T) = vecnorm(V(1:n2,I)*vv).^2*ccw'/T;
    parfor_progress(pf);
end

% Plot results
figure;
plot(1:size(B,2),B(1,:),'-k','linewidth',3)
hold on
plot(1:size(B,2),B(2,:),'--k','linewidth',2)
plot(1:size(B,2),B(3,:),':k','linewidth',2)
plot([0,size(B,2)],[1,1]*max(mu3*0.001/K),'g','linewidth',2)
legend({'$n_2=5$','$n_2=10$','$n_2=100$','$\mu_{e_1}^{(\mathrm{pp})}((-1,1))$'},'fontsize',16,'interpreter','latex','location','northeast')
xlabel('$T$','interpreter','latex','fontsize',18)
title('$\Gamma_{n_2,T}$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;




