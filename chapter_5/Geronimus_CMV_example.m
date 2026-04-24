clear
close all

% Add utils and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'utils');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'data_online');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%% CODE FOR Example in Section 5.5.1 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set parameters for the experiment

X = -pi:0.00005:pi;

a_c = @(k) 0*k+0.5i;
BETA = 2*angle(1+conj(a_c(0)));
tha = 2*asin(abs(a_c(0)));
F = @(th) (abs(th)>tha).*sqrt(cos(tha/2)^2-cos(th/2).^2)./(2*pi*abs(1+a_c(0))*abs(sin((th-BETA)/2)));
stdev = abs(abs(BETA)-abs(tha))/2;

mass = max(2/abs(1+a_c(0))^2*(abs(a_c(0)+1/2)^2-1/4),0);

%% Plot the measure

figure
plot(X,F(X),'k','linewidth',2)
hold on
plot(BETA,mass,'xk','markersize',10,'linewidth',2)
plot([BETA,BETA],[0,0.5],'--k','linewidth',1)
ax = gca; ax.FontSize = 14;
xlim([-pi,pi])
xlabel('$\theta$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

%% Form large matrix

N = 2000;
rho_c = @(k) sqrt(1-abs(a_c(k)).^2);
A = sparse(N+2,N+2);
A(1:2,1:3) = [conj(a_c(0)) conj(a_c(1))*rho_c(0) rho_c(1)*rho_c(0);
    rho_c(0) -conj(a_c(1))*a_c(0) -rho_c(1)*a_c(0)];
for j=1:round((N+2)/2)
    A([2*j+1,2*j+2],2*j:(2*j+3)) = [conj(a_c(2*j))*rho_c(2*j-1) -conj(a_c(2*j))*a_c(2*j-1) conj(a_c(2*j+1))*rho_c(2*j) rho_c(2*j)*rho_c(2*j+1);
                                    rho_c(2*j)*rho_c(2*j-1)     -rho_c(2*j)*a_c(2*j-1)     -conj(a_c(2*j+1))*a_c(2*j)  -rho_c(2*j+1)*a_c(2*j)];
end

% speed up computation for rational kernels by diagonalising polar decomposition
U = A(1:N,1:N);
[U,~,V] = svd(full(U)'); U=V*U';
[V,E] = eig(U,'vector');
c = zeros(N,1); c(1) = 1;
c2 = V\c; % change to eigenvector coordinates

%% Plot example

X = -1.5:0.0001:0;
mu = F(X);

epsilon = 0.01;

order = 2;
mu2 = 0*X;
[poles,res] = rational_kernel(order,'equi');
for mm = 1:order
    Z = exp(1i*X-1i*epsilon*poles(mm));
    mu2 = mu2 - real(res(mm)*FC(E,Z,c2))/pi/2;
end

order = 4;
mu4 = 0*X;
[poles,res] = rational_kernel(order,'equi');
for mm = 1:order
    Z = exp(1i*X-1i*epsilon*poles(mm));
    mu4 = mu4 - real(res(mm)*FC(E,Z,c2))/pi/2;
end

order = 6;
mu6 = 0*X;
[poles,res] = rational_kernel(order,'equi');
for mm = 1:order
    Z = exp(1i*X-1i*epsilon*poles(mm));
    mu6 = mu6-real(res(mm)*FC(E,Z,c2))/pi/2;
end

N = round(6/epsilon);
MU = zeros(N+1,1);
MU(1) = 1;
phi_opt4 = @(x)  1- x.^4.*(-20*abs(x).^3 + 70*x.^2 - 84*abs(x) + 35);
phi_inft = @(x) exp(-2./(1-abs(x)).*exp(-0.109550455106347./abs(x).^4));

v = c;
for j = 2:N+1
    v = U*v;
    MU(j) = v(1);
end
MU = flipud([conj(flipud(MU(2:end)));MU]/(2*pi));

nn = abs(-N:N);
g4 = phi_opt4(nn/N);
gJ = ((1-nn/N).*cos(nn/N*pi) + cot(pi/N)/N*sin(nn/N*pi));
ginf = phi_inft(nn/N);

nu2 = chebfun(gJ(:).*MU(:),[-pi pi],'trig','coeffs');
nu3 = chebfun(g4(:).*MU(:),[-pi pi],'trig','coeffs');
nu4 = chebfun(ginf(:).*MU(:),[-pi pi],'trig','coeffs');

mu(isnan(mu)) = eps;
mu(abs(mu)<10^(-14)) = eps;

figure
semilogy(X,mu,':g','linewidth',4)
hold on
semilogy(X,max(mu2,eps),'-k','linewidth',1)
xlim([-1.5,0])
ylim(10.^([-11,5]))
xlabel('$\theta$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
title('$m=2$','interpreter','latex','fontsize',30)

figure
semilogy(X,mu,':g','linewidth',4)
hold on
semilogy(X,max(mu4,eps),'-k','linewidth',1)
xlim([-1.5,0])
ylim(10.^([-11,5]))
xlabel('$\theta$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
title('$m=4$','interpreter','latex','fontsize',30)

figure
semilogy(X,mu,':g','linewidth',4)
hold on
semilogy(X,max(mu6,eps),'-k','linewidth',1)
xlim([-1.5,0])
ylim(10.^([-11,5]))
xlabel('$\theta$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
title('$m=6$','interpreter','latex','fontsize',30)

figure
semilogy(X,mu,':g','linewidth',4)
hold on
semilogy(X,max(nu2(X),eps),'-k','linewidth',1)
xlim([-1.5,0])
ylim(10.^([-11,5]))
xlabel('$\theta$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
title('Jackson','interpreter','latex','fontsize',30)

figure
semilogy(X,mu,':g','linewidth',4)
hold on
semilogy(X,max(nu3(X),eps),'-k','linewidth',1)
xlim([-1.5,0])
ylim(10.^([-11,5]))
xlabel('$\theta$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
title('Vandeven 4','interpreter','latex','fontsize',30)

figure
semilogy(X,mu,':g','linewidth',4)
hold on
semilogy(X,max(nu4(X),eps),'-k','linewidth',1)
xlim([-1.5,0])
ylim(10.^([-11,5]))
xlabel('$\theta$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
title('Bump','interpreter','latex','fontsize',30)


function d = FC(E,Z,c) % fast way to apply resolvent methods using diagonalised matrix

a = (E+Z).*conj(c);
b = (1./(E-Z)).*c;
d = sum(a.*b,1);


end







