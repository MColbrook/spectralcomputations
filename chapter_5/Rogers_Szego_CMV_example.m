close all
clear

% Add SpecSolve folder from chapter 4 (constructs the rational kernels)
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'chapter_4', 'SpecSolve');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%%% CODE FOR THE CMV MATRIX EXAMPLE %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% EXPERIMENTS IN SECTION 5.2 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set the parameters of the experiment

q = 0.5;
X = -pi:0.001:pi;
F = @(th) double(theta_func( (q),(th) )/(2*pi));
mu = F(X); % density of the measure for comparison


%% Plot cumulative distribution functions

figure % CDF for spectral measure
plot(X,cumsum(mu)/sum(mu),'k','linewidth',2)
xlim([-pi,pi])
ylim([0,1])
xlabel('$\theta$','interpreter','latex','fontsize',18)
title('Cumulative Distribution','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=22;

figure % CDF of discrete approximation
N = 20; % discretisation size
[U,~] = RS_mat(q,N);
[U,~,V]=svd(full(U)'); U=V*U';
[V,E2] = eig(U,'vector'); E3 = angle(E2(:)); c = abs(V(1,:)).^2;
[~,I] =sort(E3,'ascend');
xx = E3(I); xx = xx(:); xx = [xx(1:end)';xx(1:end)'];xx =[-pi;xx(:);pi];
yy = cumsum(c(I)'); yy = yy(:); yy = [0;yy];  yy = [yy(1:end-1)';yy(2:end)'];yy =[0;yy(:);1];
plot(xx,yy,'k','linewidth',2)
xlim([-pi,pi])
ylim([0,1])
xlabel('$\theta$','interpreter','latex','fontsize',18)
title('Discrete Approx. ($n=20$)','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=22;

figure % CDF of discrete approximation (larger N)
N = 100;
[U,~] = RS_mat(q,N);
E = eig(full(U));
[U,~,V]=svd(full(U)'); U=V*U';
[V,E2] = eig(U,'vector'); E3 = angle(E2(:)); c = abs(V(1,:)).^2;
[~,I] =sort(E3,'ascend');
xx = E3(I); xx = xx(:); xx = [xx(1:end)';xx(1:end)'];xx =[-pi;xx(:);pi];
yy = cumsum(c(I)'); yy = yy(:); yy = [0;yy];  yy = [yy(1:end-1)';yy(2:end)'];yy =[0;yy(:);1];
plot(xx,yy,'k','linewidth',2)
xlim([-pi,pi])
ylim([0,1])
xlabel('$\theta$','interpreter','latex','fontsize',18)
title('Discrete Approx. ($n=100$)','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=22;


%% Eigenvalue plots (eigenvalues computed above)

figure % without polar decomposition
plot(real(E),imag(E),'k.','markersize',16)
hold on
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k','linewidth',1)
axis equal
title('$\mathrm{Sp}(\mathcal{P}_{100}A\mathcal{P}_{100}^*)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
axis([-1.1,1.1,-1.1,1.1])

figure % with polar decomposition
plot(real(E2),imag(E2),'k.','markersize',16)
hold on
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k','linewidth',1)
axis equal
title('Sp of Unitary Part of $\mathcal{P}_{100}A\mathcal{P}_{100}^*$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
axis([-1.1,1.1,-1.1,1.1])


%%%%%%%%%%%%%%%%%%% EXPERIMENTS IN SECTION 5.4 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Errors for rational kernels

epsvec = 10.^(-2:0.1:0); % values of epsilon for test
N=200; % truncation size (in general this is larger for smaller epsilon)
[U,c] = RS_mat(q,N);

Erat = zeros(6,length(epsvec));
X = 0:0.01:pi;
mu = F(X);

for order =1:6
    [poles,res]=rational_kernel(order,'equi');
    for j=1:length(epsvec)
        epsilon=epsvec(j);
        mu1 = 0*X;
        for jj=1:length(X)
            for mm=1:order
                b = (U+exp(1i*X(jj)-1i*epsilon*poles(mm))*speye(size(U)))'*c;
                a = res(mm)*((U-exp(1i*X(jj)-1i*epsilon*poles(mm))*speye(size(U)))\c);
                mu1(jj)=mu1(jj)-real(b'*a)/pi/2;
            end
        end
        Erat(order,j) = max(mu1(:)-mu(:));
    end
end

figure
loglog(epsvec,Erat(1,:),'marker','o','MarkerSize',8,'color','k','linewidth',1.5)
hold on
loglog(epsvec,Erat(2,:),'marker','s','MarkerSize',8,'color','k','linewidth',1.5)
loglog(epsvec,Erat(3,:),'marker','^','MarkerSize',8,'color','k','linewidth',1.5)
loglog(epsvec,Erat(4,:),'marker','.','MarkerSize',25,'color','k','linewidth',1.5)
loglog(epsvec,Erat(5,:),'marker','*','MarkerSize',8,'color','k','linewidth',1.5)
loglog(epsvec,Erat(6,:),'marker','+','MarkerSize',8,'color','k','linewidth',1.5)
legend({'$m=1$','$m=2$','$m=3$','$m=4$','$m=5$','$m=6$'},'fontsize',16,'interpreter','latex','location','southeast')
xlabel('$\epsilon$','interpreter','latex','fontsize',18)
title('$\|K_\epsilon^{\bf{T}}\ast\xi_{e_1}-\rho_{e_1}\|_\infty$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;


%% Errors for trigonometric polynomial kernels

N = 200;
[U,c] = RS_mat(q,N);
MU = zeros(N+1,1);
MU(1) = 1;
epsvec = 10.^(-2:0.1:0);

phi_fejer = @(x) 1-abs(x);
phi_cosine = @(x) (1+cos(pi*x))/2;
phi_opt4 = @(x)  1- x.^4.*(-20*abs(x).^3 + 70*x.^2 - 84*abs(x) + 35);
phi_sharp_cosine = @(x) phi_cosine(x).^4.*(35-84*phi_cosine(x)+70*phi_cosine(x).^2-20*phi_cosine(x).^3);
phi_inft = @(x) exp(-2./(1-abs(x)).*exp(-0.109550455106347./abs(x).^4));

v = c;
for j=2:N+1
    v = U*v;
    MU(j) = v(1);
end
MU=[conj(flipud(MU(2:end)));MU]/(2*pi);
N0 = N;

Epoly = zeros(4,length(epsvec));
for j=1:length(epsvec)
    N=round(1/epsvec(j));
    nn = abs(-N:N);
    gF = phi_fejer(nn/N);
    gJ = (1-nn/N).*cos(nn/N*pi) + cot(pi/N)/N*sin(nn/N*pi);
    g4 = phi_opt4(nn/N);
    ginf = phi_inft(nn/N);

    nu = chebfun(gF(:).*MU((N0+1-N):(N0+1+N)),[-pi pi],'trig','coeffs');
    Epoly(1,j) = max(abs(nu(X)-mu));

    nu = chebfun(gJ(:).*MU((N0+1-N):(N0+1+N)),[-pi pi],'trig','coeffs');
    Epoly(2,j) = max(abs(nu(X)-mu));

    nu = chebfun(g4(:).*MU((N0+1-N):(N0+1+N)),[-pi pi],'trig','coeffs');
    Epoly(3,j) = max(abs(nu(X)-mu));

    nu = chebfun(ginf(:).*MU((N0+1-N):(N0+1+N)),[-pi pi],'trig','coeffs');
    Epoly(4,j) = max(abs(nu(X)-mu));
    

end

figure
loglog(round(1./epsvec),Epoly(1,:),'marker','o','MarkerSize',8,'color','k','linewidth',1.5)
hold on
loglog(round(1./epsvec),Epoly(2,:),'marker','s','MarkerSize',8,'color','k','linewidth',1.5)
loglog(round(1./epsvec),Epoly(3,:),'marker','.','MarkerSize',25,'color','k','linewidth',1.5)
loglog(round(1./epsvec),Epoly(4,:),'marker','h','MarkerSize',8,'color','k','linewidth',1.5)
legend({'Fej{\''e}r','Jackson','Vandeven 4','Bump'},'fontsize',16,'interpreter','latex','location','southwest')
xlabel('$N=\lfloor \epsilon^{-1}\rfloor$','interpreter','latex','fontsize',18)
title('$\|\xi_{N;e_1}^\sigma-\rho_{e_1}\|_\infty$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=16;




function [U,c] = RS_mat(q,N)
%% Form matrix
a_c = @(k) (-1).^k.*q.^((k+1)/2);
rho_c = @(k) sqrt(1-abs(a_c(k)).^2);
A=sparse(N+2,N+2);
A(1:2,1:3)=[conj(a_c(0)) conj(a_c(1))*rho_c(0) rho_c(1)*rho_c(0);
    rho_c(0) -conj(a_c(1))*a_c(0) -rho_c(1)*a_c(0)];
for j=1:round((N+2)/2)
    A([2*j+1,2*j+2],[2*j:2*j+3])=[conj(a_c(2*j))*rho_c(2*j-1) -conj(a_c(2*j))*a_c(2*j-1) conj(a_c(2*j+1))*rho_c(2*j) rho_c(2*j)*rho_c(2*j+1);
                                  rho_c(2*j)*rho_c(2*j-1)     -rho_c(2*j)*a_c(2*j-1)     -conj(a_c(2*j+1))*a_c(2*j)  -rho_c(2*j+1)*a_c(2*j)];
end

U=A(1:N,1:N);
c=zeros(N,1); c(1)=1;
end


function [ Y ] = theta_func( q,th )
Y=0*th;

for a=-40:40
    Y=Y+exp(-(th-2*pi*a).^2/(2*log(1/q)));
end
Y=Y*2*pi/sqrt(2*pi*log(1/q));

end









