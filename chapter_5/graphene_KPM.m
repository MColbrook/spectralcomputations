clear
close all

% Add utils and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'utils');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'data_online');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%%%%%%% CODE FOR Example 5.4.18 %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct the matrix of the operator
PHI = 1/4; % magnetic field

% construct the matrix from precomputation
load('graphene_lattice.mat')
NN = (2*n1+1)*(2*n2+1)*2;
N = 2100;

n = min(10^7,length(ord)); % number of basis sites used

d1 = spdiags(ones(2*n1,1),-1,2*n1+1,2*n1+1);
d2 = spdiags(ones(2*n2,1),-1,2*n2+1,2*n2+1);
D2 = spdiags(transpose(exp(2i*pi*PHI*(-n1:1:n1))),0,2*n1+1,2*n1+1);

T1 = [sparse(round(NN/2),round(NN/2)), speye(round(NN/2));
      speye(round(NN/2)), sparse(round(NN/2),round(NN/2))];
T2 = [sparse(round(NN/2),round(NN/2)), kron(speye(2*n2+1),d1');
      kron(speye(2*n2+1),d1),sparse(round(NN/2),round(NN/2)) ];
T3 = [sparse(round(NN/2),round(NN/2)), kron(d2',D2);
      kron(d2,D2'),sparse(round(NN/2),round(NN/2)) ];

H0 = T1+T2+T3;
H0 = H0(ord,ord); % reorder for operator on l^2(N)

%% Rescale operator to have spectrum in [-1,1]

SC = eigs(H0(1:10^4,1:10^4),1,2.6)*1.00001; % scaling to make spectrum in (-1,1)
H0 = H0/SC;
SC = eigs(H0(1:10^6,1:10^6),1,1);
H0 = H0/SC;

%% Compute the moments (uncomment to perform computation or use saved data)

% MU = zeros(N+1,1);
% MU(1) = 1;
% H0 = H0(1:n,1:n);
% 
% b=sparse(n,1);           b(1)=1;
% v0 = b;
% v1 = H0*b;
% MU(2) = v1(1);
% for j=3:N+1
%     a = 2*H0*v1-v0;
%     v0 = v1;
%     v1 = a;
%     MU(j) = a(1);
%     if j==1001
%         sum(abs(a)>0)
%     end
%     if j==2001
%         sum(abs(a)>0)
%     end
% end
% MU(2:end) = 2*MU(2:end);

load('graph_moments.mat')

%% Compute measure for comparison

tt = length(MU)-1;
phi_opt4 = @(x)  1- x.^4.*(-20*abs(x).^3 + 70*x.^2 - 84*abs(x) + 35);
phi_inft = @(x) exp(-2./(1-abs(x)).*exp(-0.109550455106347./abs(x).^4));
nn = 0:tt;
gJ = ((1-nn/tt).*cos(nn/tt*pi) + cot(pi/tt)/tt*sin(nn/tt*pi)); gJ = gJ(:);
ginf = phi_inft(nn/tt);
ginf = ginf(:);
nu0 = chebfun(real(ginf(1:(tt+1)).*MU(1:(tt+1))),'coeffs'); % reference measure

%% Compute the errors in N

X0 = 0; X1 = 0.12; % points where we compute error
mu0 = 0; % actual value
mu1 = nu0(X1);
Nvec = 100:2020;

E0 = zeros(4,length(Nvec)); E1 = E0;
for jj=1:length(Nvec)
    nn = 0:Nvec(jj);
    g4 = phi_opt4(nn/Nvec(jj));
    gJ = ((1-nn/Nvec(jj)).*cos(nn/Nvec(jj)*pi) + cot(pi/Nvec(jj))/Nvec(jj)*sin(nn/Nvec(jj)*pi));
    ginf = phi_inft(nn/Nvec(jj));
    gJ = gJ(:);
    g4 = g4(:);

    ginf = ginf(:);
    
    nu1 = chebfun(real(MU(1:(Nvec(jj)+1))),'coeffs');
    nu2 = chebfun(real(gJ(1:(Nvec(jj)+1)).*MU(1:(Nvec(jj)+1))),'coeffs');
    nu3 = chebfun(real(g4(1:(Nvec(jj)+1)).*MU(1:(Nvec(jj)+1))),'coeffs');
    nu4 = chebfun(real(ginf(1:(Nvec(jj)+1)).*MU(1:(Nvec(jj)+1))),'coeffs');
    E0(1,jj) = abs(nu1(X0)-mu0);
    E0(2,jj) = abs(nu2(X0)-mu0);
    E0(3,jj) = abs(nu3(X0)-mu0);
    E0(4,jj) = abs(nu4(X0)-mu0);

    E1(1,jj) = abs(nu1(X1)-mu1);
    E1(2,jj) = abs(nu2(X1)-mu1);
    E1(3,jj) = abs(nu3(X1)-mu1);
    E1(4,jj) = abs(nu4(X1)-mu1);

end
E0 = E0/pi;
E1 = E1/pi;

% Plot the results
figure
loglog(Nvec,E0(1,:),'-k','linewidth',1)
hold on
loglog(Nvec,E0(2,:),'-k','linewidth',3)
loglog(Nvec,E0(3,:),'--k','linewidth',1)
loglog(Nvec,E0(4,:),':k','linewidth',1)
legend({'No kernel','Jackson','Vandeven 4','Bump'},'fontsize',16,'interpreter','latex','location','northwest')
ax=gca; ax.FontSize=16;
xlabel('$N$','interpreter','latex','fontsize',18)
title('$E_{N,\infty}^{\rm KPM}, x_0=0$','interpreter','latex','fontsize',18)

figure
loglog(Nvec,E1(1,:),'-k','linewidth',1)
hold on
loglog(Nvec,E1(2,:),'-k','linewidth',3)
loglog(Nvec,E1(3,:),'--k','linewidth',1)
loglog(Nvec,E1(4,:),':k','linewidth',1)
ax=gca; ax.FontSize=16;
xlabel('$N$','interpreter','latex','fontsize',18)
title('$E_{N,\infty}^{\rm KPM},x_0=0.12$','interpreter','latex','fontsize',18)


%% Plot the approximations at N = 1000

N = 1000;
nn = 0:N;
g4 = phi_opt4(nn/N);
gJ = ((1-nn/N).*cos(nn/N*pi) + cot(pi/N)/N*sin(nn/N*pi));
ginf = phi_inft(nn/N);
gJ = gJ(:);
g4 = g4(:);

ginf = ginf(:);

nu1 = chebfun(real(MU(1:(N+1))),'coeffs');
nu2 = chebfun(real(gJ(1:(N+1)).*MU(1:(N+1))),'coeffs');
nu3 = chebfun(real(g4(1:(N+1)).*MU(1:(N+1))),'coeffs');
nu4 = chebfun(real(ginf(1:(N+1)).*MU(1:(N+1))),'coeffs');

f=figure;
tiledlayout(2,2,'TileSpacing','Compact');
nexttile
dummyh = plot(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
hold on
plot(nu1,'-k');
title('No kernel','fontsize',16,'interpreter','latex')
ylim([-2,15])
xlim([0,1.05])
xlabel('$x$','interpreter','latex','fontsize',16)
grid minor
nexttile
plot(max(nu2,0),'-k')
title('Jackson','fontsize',16,'interpreter','latex')
ylim([-2,15])
xlim([0,1.05])
xlabel('$x$','interpreter','latex','fontsize',16)
grid minor
nexttile
plot(max(0,nu3),'-k')
title('Vandeven 4','fontsize',16,'interpreter','latex')
ylim([-2,15])
xlim([0,1.05])
xlabel('$x$','interpreter','latex','fontsize',16)
grid minor
nexttile
plot(max(0,nu4),'-k')
title('Bump','fontsize',16,'interpreter','latex')
ylim([-2,15])
xlim([0,1.05])
xlabel('$x$','interpreter','latex','fontsize',16)
grid minor
f.Position = 1.0e+03*[0.1000    0.0977    1.1200    0.5160];


