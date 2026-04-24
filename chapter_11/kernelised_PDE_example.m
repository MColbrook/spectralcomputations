clear
close all

% Add Koopman algorithms and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'Koopman_algorithms');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'Koopman_datasets');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%% CODE FOR FITZHUGH-NAGUMO PDE EXAMPLE %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1) % set random seed for reproducibility

%% Set the parameters
% grid for pseudospectra
x_pts=0:0.005:1.2;    y_pts=-0.01:0.005:0.8;
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);

N = 100; % size of dictionary used

%% Create the data set
nn = 64;
sigma = 0.1;
L = 20;
ns = 50;

Xpts = (((1:nn)-1/2)/nn*L); Xpts = Xpts(:);
D = -((0:nn-1)*pi/L).^2; D = D(:);
E = cos(Xpts*(0:nn-1)*pi/L);
EM = inv(E);
delta_t = 1;

ODEFUN=@(t,y) [D.*y(1:nn)+y(1:nn)-y((nn+1):2*nn)-EM*((E*y(1:nn)).^3);
                    4*D.*y((nn+1):2*nn)+0.02*(y(1:nn)-2*y((nn+1):2*nn)+0.03)];
options = odeset('RelTol',1e-4,'AbsTol',1e-4);

U0=[-0.12+zeros(nn/2,1);0.12+zeros(nn/2,1);-0.02+zeros(nn/2,1);0.02+zeros(nn/2,1)];
U0=[EM zeros(nn,nn);zeros(nn,nn) EM]*U0;
[~,U]=ode45(ODEFUN,[0.000001 3000],U0,options);
U = U';
U0 = U(:,end);

X = []; Y = [];

for j = 1:100
    [~,U]=ode45(ODEFUN,[0.000001 (delta_t*(1:ns))],U0,options);
    U = U';
    U = [E zeros(nn,nn);zeros(nn,nn) E]*U;
    X = [X,U0,U(:,2:end-1)];
    Y = [Y,U(:,2:end)];

    u1 = randn(1)*sigma; u2 = randn(1)*sigma; u3 = randn(1)*sigma;
    U0 = U(:,end);
    U0(1:nn) = U0(1:nn) + u1*exp(-(Xpts-7.5).^2) + u2*exp(-(Xpts-10).^2) + u3*exp(-(Xpts-12.5).^2);
    U0=[EM zeros(nn,nn);zeros(nn,nn) EM]*U0;
end

%% Run kEDMD and kResDMD

% Polynomial kernel 1
[G,K,L] = kResDMD(X,Y,'N',N,'type',"Linear");
[V,LAM] = eig(K,'vector');
yp1 = real(sqrt(dot(V,L*V+V*diag(abs(LAM)).^2-K'*V*diag(LAM)-K*V*diag(conj(LAM)))./dot(V,V))); % dual residual
[~,I]= sort(abs(1-abs(LAM)),'ascend');
xp1 = LAM(I); yp1 = yp1(I);
RESp1 = KoopPseudoSpec(G,K,L,z_pts,'Parallel','off');	% compute pseudospectra
RESp1 = reshape(RESp1,length(y_pts),length(x_pts));

% Polynomial kernel 50
[G,K,L] = kResDMD(X,Y,'N',N,'type',50);
[V,LAM] = eig(K,'vector');
yp50 = real(sqrt(dot(V,L*V+V*diag(abs(LAM)).^2-K'*V*diag(LAM)-K*V*diag(conj(LAM)))./dot(V,V))); % dual residual
[~,I]= sort(abs(1-abs(LAM)),'ascend');
xp50 = LAM(I); yp50 = yp50(I);
RESp50 = KoopPseudoSpec(G,K,L,z_pts,'Parallel','off');	% compute pseudospectra
RESp50 = reshape(RESp50,length(y_pts),length(x_pts));

% Gauss kernel
[G,K,L] = kResDMD(X,Y,'N',N);
[V,LAM] = eig(K,'vector');
yg = real(sqrt(dot(V,L*V+V*diag(abs(LAM)).^2-K'*V*diag(LAM)-K*V*diag(conj(LAM)))./dot(V,V))); % dual residual
[~,I]= sort(abs(1-abs(LAM)),'ascend');
xg = LAM(I); yg = yg(I);
RESg = KoopPseudoSpec(G,K,L,z_pts,'Parallel','off');	% compute pseudospectra
RESg = reshape(RESg,length(y_pts),length(x_pts));

% Sobolev kernel
d = mean(vecnorm(X-mean(X,2)));
nr = @(x,y) sqrt(-2*real(y'*x)+dot(x,x)+dot(y,y)')/d;
kernel_f = @(x,y) (1+nr(x,y)).*exp(-nr(x,y));
[G,K,L] = kResDMD(X,Y,'N',N,'kernel',kernel_f);
[V,LAM] = eig(K,'vector');
ys = real(sqrt(dot(V,L*V+V*diag(abs(LAM)).^2-K'*V*diag(LAM)-K*V*diag(conj(LAM)))./dot(V,V))); % dual residual
[~,I]= sort(abs(1-abs(LAM)),'ascend');
xs = LAM(I); ys = ys(I);
RESs = KoopPseudoSpec(G,K,L,z_pts,'Parallel','off');	% compute pseudospectra
RESs = reshape(RESs,length(y_pts),length(x_pts));

%% Plot the eigenvalues

figure
tiledlayout(2,1,'TileSpacing','tight')
nexttile
plot(angle(xp1)/(delta_t),log(abs(xp1))/(delta_t),'k.','markersize',18)
ylim([-0.03,0.001])
xlim([-0.5,0.5])
title('Linear','interpreter','latex','fontsize',18)
ylabel('$\log(|\lambda|)/\Delta t$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

nexttile
plot(angle(xp50)/(delta_t),log(abs(xp50))/(delta_t),'k.','markersize',18)
ylim([-0.03,0.001])
xlim([-0.5,0.5])
title('Polynomial, $n=50$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{arg}(\lambda)/\Delta t$','interpreter','latex','fontsize',18)
ylabel('$\log(|\lambda|)/\Delta t$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;


figure
tiledlayout(2,1,'TileSpacing','tight')
nexttile
plot(angle(xg)/(delta_t),log(abs(xg))/(delta_t),'k.','markersize',18)
ylim([-0.03,0.001])
xlim([-0.5,0.5])
title('Gaussian','interpreter','latex','fontsize',18)
ylabel('$\log(|\lambda|)/\Delta t$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

nexttile
plot(angle(xs)/(delta_t),log(abs(xs))/(delta_t),'k.','markersize',18)
ylim([-0.03,0.001])
xlim([-0.5,0.5])
title('Mat\''ern, $\nu = 3/2$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{arg}(\lambda)/\Delta t$','interpreter','latex','fontsize',18)
ylabel('$\log(|\lambda|)/\Delta t$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

%% Plot pseudospectra

figure
hold on
v=(10.^(-20:0.2:0));
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RESp1)),log10(v));
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),-reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RESp1)),log10(v));
cbh=colorbar;
cbh.Ticks=log10(10.^(-10:1:0));
cbh.TickLabels=["1e-10","1e-9","1e-8","1e-7","1e-6","1e-5","1e-5","1e-3","1e-2","1e-1","1"];
clim([-12,0]);
reset(gcf)
set(gca,'YDir','normal')
colormap gray
title('Linear','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18; axis equal;  axis([0,1.2,-0.8,0.8])
hold on
plot(real(xp1),imag(xp1),'.k','markersize',12);
box on

figure
hold on
v=(10.^(-20:0.25:0));
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RESp50)),log10(v));
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),-reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RESp50)),log10(v));
cbh=colorbar;
cbh.Ticks=log10(10.^(-10:1:0));
cbh.TickLabels=["1e-10","1e-9","1e-8","1e-7","1e-6","1e-5","1e-5","1e-3","1e-2","1e-1","1"];
clim([-6,0]);
reset(gcf)
set(gca,'YDir','normal')
colormap gray
title('Polynomial, $n=50$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18; axis equal;  axis([0,1.2,-0.8,0.8])
hold on
plot(real(xp50),imag(xp50),'.k','markersize',12);
box on

figure
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RESg)),log10(v));
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),-reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RESg)),log10(v));
cbh=colorbar;
cbh.Ticks=log10(10.^(-10:1:0));
cbh.TickLabels=["1e-10","1e-9","1e-8","1e-7","1e-6","1e-5","1e-5","1e-3","1e-2","1e-1","1"];
clim([-6,0]);
reset(gcf)
set(gca,'YDir','normal')
colormap gray
title('Gaussian','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18; axis equal;  axis([0,1.2,-0.8,0.8])
hold on
plot(real(xg),imag(xg),'.k','markersize',12);
box on

figure
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RESs)),log10(v));
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),-reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RESs)),log10(v));
cbh=colorbar;
cbh.Ticks=log10(10.^(-10:1:0));
cbh.TickLabels=["1e-10","1e-9","1e-8","1e-7","1e-6","1e-5","1e-5","1e-3","1e-2","1e-1","1"];
clim([-6,0]);
reset(gcf)
set(gca,'YDir','normal')
colormap gray
title('Mat\''ern, $\nu = 3/2$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18; axis equal;  axis([0,1.2,-0.8,0.8])
hold on
plot(real(xs),imag(xs),'.k','markersize',12);
box on



