clear
close all

% Add Koopman algorithms and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'Koopman_algorithms');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'Koopman_datasets');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%%% CODE FOR DAMPED DUFFING EXAMPLE %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1) % set random seed for reproducibility

%% Set parameters
M1 = 10^4; M2 = 5;
M = M1*M2; % number of data points
delta_t = 0.3; % time step
ODEFUN = @(t,y) [y(2);-0.3*y(2)+y(1)-y(1).^3];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
N = 500;
PHI = @(r) exp(-r*2); % radial basis function used (others also work well)

%% Produce the centres
X = zeros(2,M);
Y=[];
for jj=1:M1
    Y0=(rand(2,1)-0.5)*4;
    [~,Y1]=ode45(ODEFUN,[0 0.000001 (1:(3+M2))*delta_t],Y0,options);
    Y1=Y1';
    X(:,(jj-1)*M2+1:jj*M2) = Y1(:,[1,3:M2+1]);
    Y(:,(jj-1)*M2+1:jj*M2) = Y1(:,3:M2+2);
end
[~,C] = kmeans([X';Y'],N,'MaxIter',500); % find centers

%% Produce the snapshots
X = zeros(2,M);
Y=[];
for jj=1:M1
    Y0=(rand(2,1)-0.5)*4;
    [~,Y1]=ode45(ODEFUN,[0 0.000001 (1:(3+M2))*delta_t],Y0,options);
    Y1=Y1';
    X(:,(jj-1)*M2+1:jj*M2) = Y1(:,[1,3:M2+1]);
    Y(:,(jj-1)*M2+1:jj*M2) = Y1(:,3:M2+2);
end

d=mean(vecnorm(X-mean(X,2))); % scaling for radial function
PX = zeros(M,N); PY = zeros(M,N);

for j = 1:N
    R = sqrt((X(1,:)-C(j,1)).^2+(X(2,:)-C(j,2)).^2);
    PX(:,j) = PHI(R(:)/d);
    R = sqrt((Y(1,:)-C(j,1)).^2+(Y(2,:)-C(j,2)).^2);
    PY(:,j) = PHI(R(:)/d);
end

%% EDMD
K = PX(1:M,:)\PY(1:M,:);
[V,LAM] = eig(K,'vector');

figure % plot EDMD eigenvalues
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-k')
hold on
plot(real(LAM),imag(LAM),'k.','markersize',18)
[~,I]= sort(abs((0.897+0.423i)-LAM),'ascend');
plot(real(LAM(I(1))),imag(LAM(I(1))),'ko','markersize',12)
plot(real(LAM(I(1))),-imag(LAM(I(1))),'ko','markersize',12)
[~,I]= sort(abs((0.667+0.722i)-LAM),'ascend');
plot(real(LAM(I(1))),imag(LAM(I(1))),'ko','markersize',12)
plot(real(LAM(I(1))),-imag(LAM(I(1))),'ko','markersize',12)
[~,I]= sort(abs((0.378+0.867i)-LAM),'ascend');
plot(real(LAM(I(1))),imag(LAM(I(1))),'ko','markersize',12)
plot(real(LAM(I(1))),-imag(LAM(I(1))),'ko','markersize',12)
axis equal
axis([-1.15,1.15,-1.15,1.15])
title(sprintf('Eigenvalues ($N=%d$, damped)',N),'interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

return
%% plot some of the eigenfunctions corresponding to lambda = 1 (uncomment each line for different eigenfunctions)
% [~,I]= sort(abs((1)-LAM),'ascend');
% [~,I]= sort(abs((0.897+0.423i)-LAM),'ascend');
% [~,I]= sort(abs((0.667+0.722i)-LAM),'ascend');
[~,I]= sort(abs((0.378+0.867i)-LAM),'ascend');

LAM = LAM(I); V = V(:,I);
u=real([PX;PY]*V(:,1)); u = u - mean(u);
figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter([X(1,:),Y(1,:)],[X(2,:),Y(2,:)],5,u,'filled');
tt=[X(2,:),Y(2,:)];
c=colormap(bone);
clim([-2*std(u)+mean(u),2*std(u)+mean(u)])
colorbar('south')
axis(axes1,'tight'); grid off; axis off
hold on
ylim([min(tt(:))-1,max(tt(:))])
plot(1,0,'wx','markersize',20,'linewidth',2)
set(axes1,'DataAspectRatio',[1 1 1]);

%% ResDMD algorithms
R = vecnorm(PY*V-PX*V*diag(LAM))./vecnorm(PX*V); % residuals of the EDMD eigenpairs, book uses N =250 (adjust above)

figure % plot residuals of EDMD eigenvalues
loglog([min(R),1],[min(R),1],':k','linewidth',1)
loglog([0.001,1],[0.001,1],':k','linewidth',1)
hold on
loglog(sqrt(abs(abs(LAM).^2-1)),R,'k.','markersize',20)
xlabel('$\sqrt{1-|\lambda|^2}$','interpreter','latex','fontsize',18)
title('EDMD Eigenpair Residuals (damped)','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
xlim([0.01,1])

% compute pseudospectra
x_pts = -1.5:0.04:1.5;    y_pts = -0.04:0.04:1.5; % use conjuaget symmetry to restruct y values
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);
RES = KoopPseudoSpecQR(PX,PY,1/M,z_pts,'Parallel','off'); RES = reshape(RES,length(y_pts),length(x_pts));

% plot pseudospectra
v=(10.^(-4:0.1:1));
figure
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(max(min(v),real(RES))),log10(v));
hold on % use conjugate symmetry to reduce number of grid points to test
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),-reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(max(min(v),real(RES))),log10(v));
cbh=colorbar;
cbh.Ticks=log10([0.005,0.01,0.1,1]);
cbh.TickLabels=[0,0.01,0.1,1];
clim([-2,0]);
reset(gcf)
set(gca,'YDir','normal')
colormap gray
hold on
plot(real(LAM),imag(LAM),'k.','markersize',12)
axis equal tight;  axis([x_pts(1),x_pts(end),-y_pts(end),y_pts(end)])
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
title(sprintf('$\\mathrm{Sp}_\\epsilon(\\mathcal{K})$ ($N=%d$, damped)',N),'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

