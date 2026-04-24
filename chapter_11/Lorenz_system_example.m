clear
close all

% Add Koopman algorithms and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'Koopman_algorithms');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'Koopman_datasets');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%%% CODE FOR LORENZ SYSTEM EXAMPLE %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1) % set random seed for reproducibility

%% Set parameters for numerical integration of ODE
options = odeset('RelTol',1e-14,'AbsTol',1e-14);
SIGMA = 10;   BETA = 8/3;   RHO = 28;
ODEFUN = @(t,y) [SIGMA*(y(2)-y(1));y(1).*(RHO-y(3))-y(2);y(1).*y(2)-BETA*y(3)];

%% Set parameters for the example
N = 1200;    % number of delays (this is an upper bound, below we reduce at various stages)
M = 10^4;    % number of snapshots      
dt = 0.05;   % time step for trajectory sampling

%% Produce the data
M2 = 10^5; % for visualisation of functions we use more data points
Y0 = (rand(3,1)-0.5)*4; % initial point off the Lorenz attractor
[~,Y0] = ode45(ODEFUN,[0.000001 1, 100],Y0,options); Y0 = Y0(end,:)'; % new initial point on the Lorenz attractor
h = 1; % number of time steps for delay
[~,DATA] = ode45(ODEFUN,[0.000001 (1:((M2+h*(N+1))))*dt],Y0,options);

%% Use delay embedding as the dictionary
PX1=zeros(M2,N); PX1(:,1)=DATA(1:M2,1);
PX2=zeros(M2,N); PX2(:,1)=DATA(1:M2,2);
PX3=zeros(M2,N); PX3(:,1)=DATA(1:M2,3);
PY1=zeros(M2,N); PY1(:,1)=DATA((1:M2)+1,1);
PY2=zeros(M2,N); PY2(:,1)=DATA((1:M2)+1,2);
PY3=zeros(M2,N); PY3(:,1)=DATA((1:M2)+1,3);

for j=2:N
    PX1(:,j)=DATA((1:M2)+h*(j-1),1);
    PX2(:,j)=DATA((1:M2)+h*(j-1),2);
    PX3(:,j)=DATA((1:M2)+h*(j-1),3);
    PY1(:,j)=DATA((1:M2)+1+h*(j-1),1);
    PY2(:,j)=DATA((1:M2)+1+h*(j-1),2);
    PY3(:,j)=DATA((1:M2)+1+h*(j-1),3);
end

%% Plot EDMD eigenvalues and eigenfunctions
N = 100; % 25, 50, 100 in book

PX = [PX1(1:M,1:N),PX2(1:M,1:N),PX3(1:M,1:N)];
PY = [PY1(1:M,1:N),PY2(1:M,1:N),PY3(1:M,1:N)];

K = PX\PY;
[V,LAM] = eig(K,'vector');
[~,I]=sort(abs(LAM-1),'ascend');   
V =V(:,I); LAM =LAM(I);
R = vecnorm(PY*V-PX*V*diag(LAM))./vecnorm(PX*V); % residuals of the EDMD eigenpairs

figure % plot EDMD eigenvalues
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-k')
hold on
plot(real(LAM),imag(LAM),'k.','markersize',18)
axis equal
axis([-1.15,1.15,-1.15,1.15])
title(sprintf('EDMD Eigenvalues ($N=%d$)',N),'interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

LL = 8.25; % spectral parameter for plot
w = log(LAM)/(1i*dt);
I = find(abs(w-LL)==min(abs(w-LL)));

figure % plot EDMD eigenfunction
u = [PX1(:,1:N),PX2(:,1:N),PX3(:,1:N)]*V([(1:N),(1:N)+N,(1:N)+2*N],I(1));
u = u(1:min(10^6,M2),:);
u = real(u/mean(u(1:300)));
u = real(u/mean(u(1:300)));
[~,II]=sort(u,'ascend');
scatter3(PX1(II,1),PX2(II,1),PX3(II,1),6,u(II),'filled');
c=colormap(gray);
c = brighten(c(1:end,:),0.4); colormap(c)
clim([-2*std(u)+mean(u),2*std(u)+mean(u)])
view(gca,[13.1786087602293 -1.28469255513244]);
axis tight; grid off; axis off
axis equal tight
set(gca,'DataAspectRatio',[1 1 1]);
xlabel('$X$','interpreter','latex','fontsize',14)
ylabel('$Y$','interpreter','latex','fontsize',14)
zlabel('$Z$','interpreter','latex','fontsize',14,'rotation',0)
t1 = sprintf('EDMD, $N=%d$,',N);
t2 = sprintf('$\\lambda\\approx e^{\\mathrm{i}8.25\\Delta t}$, $\\mathrm{res}=%.3f$',R(I(1)));
title({t1,t2},'interpreter','latex','fontsize',18)

%% ResDMD algorithms
R = vecnorm(PY*V-PX*V*diag(LAM))./vecnorm(PX*V); % residuals of the EDMD eigenpairs

figure % plot residuals of EDMD eigenvalues
loglog([0.006,1],[0.006,1],':k','linewidth',1)
hold on
loglog(sqrt(abs(abs(LAM).^2-1)),R,'k.','markersize',20)
xlabel('$\sqrt{1-|\lambda|^2}$','interpreter','latex','fontsize',18)
title('EDMD Eigenpair Residuals','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

LL = 1;
[~,R,V2] = KoopPseudoSpecQR(PX,PY,1/M,[],'z_pts2',exp(1i*LL*dt)); % compute pseudoeigenfunction

figure % plot pseudoeigenfunction
u = [PX1(:,1:N),PX2(:,1:N),PX3(:,1:N)]*V2;
u = u(1:min(10^6,M2),:);
u = real(u/mean(u(1:300)));
u = real(u/mean(u(1:300)));
[~,I]=sort(u,'ascend');
scatter3(PX1(I,1),PX2(I,1),PX3(I,1),6,u(I),'filled');
colormap(brighten(brewermap([],'RdYlBu'),-0.3))
c=colormap(gray);
c = brighten(c(1:end,:),0.4); colormap(c)
clim([-2*std(u)+mean(u),2*std(u)+mean(u)])
view(gca,[13.1786087602293 -1.28469255513244]);
axis tight; grid off; axis off
axis equal tight
set(gca,'DataAspectRatio',[1 1 1]);
xlabel('$X$','interpreter','latex','fontsize',14)
ylabel('$Y$','interpreter','latex','fontsize',14)
zlabel('$Z$','interpreter','latex','fontsize',14,'rotation',0)
t1 = sprintf('ResDMD, $N=%d$,',N);
t2 = sprintf('$\\lambda\\approx e^{\\mathrm{i}\\Delta t}$, $\\mathrm{res}=%.3f$',R);
title({t1,t2},'interpreter','latex','fontsize',18)

LL = 8.25;
[~,R,V2] = KoopPseudoSpecQR(PX,PY,1/M,[],'z_pts2',exp(1i*LL*dt)); % compute pseudoeigenfunction

figure % plot pseudoeigenfunction
u = [PX1(:,1:N),PX2(:,1:N),PX3(:,1:N)]*V2;
u = u(1:min(10^6,M2),:);
u = real(u/mean(u(1:300)));
u = real(u/mean(u(1:300)));
[~,I]=sort(u,'ascend');
scatter3(PX1(I,1),PX2(I,1),PX3(I,1),6,u(I),'filled');
colormap(brighten(brewermap([],'RdYlBu'),-0.3))
c=colormap(gray);
c = brighten(c(1:end,:),0.4); colormap(c)
clim([-2*std(u)+mean(u),2*std(u)+mean(u)])
view(gca,[13.1786087602293 -1.28469255513244]);
axis tight; grid off; axis off
axis equal tight
set(gca,'DataAspectRatio',[1 1 1]);
xlabel('$X$','interpreter','latex','fontsize',14)
ylabel('$Y$','interpreter','latex','fontsize',14)
zlabel('$Z$','interpreter','latex','fontsize',14,'rotation',0)
t1 = sprintf('ResDMD, $N=%d$,',N);
t2 = sprintf('$\\lambda\\approx e^{\\mathrm{i}8.25\\Delta t}$, $\\mathrm{res}=%.3f$',R);
title({t1,t2},'interpreter','latex','fontsize',18)

%% mpEDMD algorithm

N = 50;
g = @(x,y,z) tanh((x.*y-3*z)/5);    % observable

PX = g(PX1(1:M,1:N),PX2(1:M,1:N),PX3(1:M,1:N));
PY = g(PY1(1:M,1:N),PY2(1:M,1:N),PY3(1:M,1:N));
[~,mpV,mpD] = mpEDMDqr(PX,PY,1/M);
LAM = diag(mpD);

figure % plot mpEDMD eigenvalues
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-k')
hold on
plot(real(LAM),imag(LAM),'k.','markersize',18)
axis equal
axis([-1.15,1.15,-1.15,1.15])
title(sprintf('mpEDMD Eigenvalues ($N=%d$)',N),'interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

% Compute spectral measure
c = zeros(N,1); c(1) = 1; % coefficients of g in Krylov basis
piE = diag(mpD); TH=angle(piE*exp(1i*eps)); % eps here is to take into account MATLAB's convention for angle
G = (PX'*PX)/M;
MU = abs(mpV'*G*c).^2;

% Cdf plots
[~,Ib] = sort(TH(:),'ascend');
THp=TH(Ib); THp=[THp(:)-10^(-14),THp(:)]'; THp=THp(:);
cdf=0*THp;
cc=0;
for j=1:length(TH)
    cdf(2*j-1)=cc;    cc=cc+MU(Ib(j));    cdf(2*j)=cc;
end

THp = [-pi;THp(:);pi]; cdf = [0;cdf(:);sum(MU)]; % for visualisation

figure
plot(THp,cdf/sum(MU),'k','linewidth',2)
ylim([0,1]); xlim([-pi,pi]);
title(sprintf('mpEDMD CDF ($N=%d$)',N),'interpreter','latex','fontsize',18)
xlabel('$\theta$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

%% mpEDMD convergence

% reference measure
N = 150;
g = @(x,y,z) tanh((x.*y-3*z)/5);    % observable

PX = g(PX1(1:M,1:N),PX2(1:M,1:N),PX3(1:M,1:N));
PY = g(PY1(1:M,1:N),PY2(1:M,1:N),PY3(1:M,1:N));
[~,mpV,mpD] = mpEDMDqr(PX,PY,1/M);
G = (PX'*PX)/M;

c = zeros(N,1); c(1) = 1; % coefficients of g in Krylov basis
TH1=angle(diag(mpD));
MU1=abs(mpV'*G*c).^2;
MU1 = MU1/sum(MU1);

% compute errors for different N
nvec = 1:1:100;
Z = zeros(length(nvec),1);

for jj = 1:length(nvec)
    n = nvec(jj)
    [~,mpV,mpD] = mpEDMDqr(PX(:,1:n),PY(:,1:n),1/M);
    c = zeros(n,1); c(1) = 1; % coefficients of g in Krylov basis
    TH2=angle(diag(mpD));
    MU2=abs(mpV'*G(1:n,1:n)*c).^2;
    MU2 = MU2/sum(MU2);
    Z(jj) = W1(TH1,MU1,TH2,MU2,-4:0.00005:4);
end

% cplot results
figure
loglog(nvec,Z,'k','linewidth',2)
hold on
loglog(nvec,pi./(nvec),'--k','linewidth',1)

ylim([0.01,1])
axis tight

text(10,0.6,'$\pi/N$ bound','interpreter','latex','fontsize',18,'rotation',-30);
title('$W_1(\xi_{g},\xi_{\mathbf{g}}^{(N,M)})$','interpreter','latex','fontsize',18)
xlabel('$N$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;





function [Z] = W1(TH1,MU1,TH2,MU2,alpha)
% Computes Wasserstein-1 distance of two discrete measures on periodic interval

break_pts=[TH1(:);TH2(:)];
ind=[TH1(:)*0;TH2(:)*0+1];
MU=[MU1(:);MU2(:)];

[break_pts,I] = sort(break_pts,'ascend');
ind=ind(I); MU=MU(I);

L=length(I);

F1 = zeros(L,1); F2 = F1;

if ind(1)==0
    F1(1) = MU(1);
else
    F2(1) = MU(1);
end

for j=2:L
    
    if ind(j)==0
        F1(j) = F1(j-1) + MU(j);
        F2(j) = F2(j-1);
    else
        F2(j) = F2(j-1) + MU(j);
        F1(j) = F1(j-1);
    end
end
gaps = break_pts(2:L)-break_pts(1:L-1);
gaps = transpose(gaps(:));


Z = 10;
alpha = transpose(alpha(:));
Z = abs(alpha)*(break_pts(1)+pi+pi-break_pts(end)) + ...
    gaps*abs(F1(1:L-1)-F2(1:L-1)-alpha);
Z = min(Z);

end
