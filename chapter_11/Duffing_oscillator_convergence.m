clear
close all

% Add Koopman algorithms and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'Koopman_algorithms');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'Koopman_datasets');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%%%%%%%% CODE FOR FIGURE 11.6 %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1) % set random seed for reproducibility

%% Set parameters
M1 = 10^4; M2 = 5;
M = M1*M2; % number of data points
delta_t = 0.3; % time step
ODEFUN = @(t,y) [y(2);y(1)-y(1).^3];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
N = 1000;
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
p = randperm(N); C = C(p,:);

%% Produce the snapshots
M1 = 3*M1; M = M1*M2;
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

%% Convergence in N
Nvec = round(10.^(0:0.05:3));
g = @(x,y) exp(-(x+y).^2);
gg = g(X(1,:)',X(2,:)');
gg2 = g(Y(1,:)',Y(2,:)'); % one time step forward
errN = zeros(length(Nvec),1);

for jj = 1:length(Nvec)
    N = Nvec(jj)
    xi = PX(:,1:N)\gg;
    K = PX(:,1:N)\PY(:,1:N);
    errN(jj) = norm(PX(1:5:end,1:N)*K*xi-gg2(1:5:end)); % 1:5 to avoid the trajectories and sample usual Lebesgue measure
end

figure
loglog(Nvec,errN/norm(gg(1:5:end)),'k','linewidth',2)
xlabel('$N$','interpreter','latex','fontsize',18)
title('$\|\mathcal{K}g-\mathcal{P}_N^*\mathcal{P}_N\mathcal{K}g_N\|/\|g\|$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

%% Convergence in M
Mvec = round(10.^(3:0.02:5)/5)*5;
Mvec = Mvec(Mvec<M);
errM = zeros(length(Mvec),3);

ct = 1;
for NNN = [100,250,500]
    Kinf = PX(1:end,1:NNN)\PY(1:end,1:NNN);
    for jj = 1:length(Mvec)
        Mvec(jj)
        K = PX(1:Mvec(jj),1:NNN)\PY(1:Mvec(jj),1:NNN);
        errM(jj,ct) = norm(K-Kinf,'fro')./norm(Kinf,'fro');
    end
    ct = ct + 1;
end

figure
loglog(Mvec,errM(:,1),'k','linewidth',2)
hold on
loglog(Mvec,errM(:,2),'--k','linewidth',2)
loglog(Mvec,errM(:,3),':k','linewidth',2)
xlabel('$M$','interpreter','latex','fontsize',18)
ylim([0.01,1])
legend({'$N=100$','$N=250$','$N=500$'},'fontsize',16,'interpreter','latex','location','northeast')
title('$\|\mathbf{K}_{N,M}-\mathbf{K}_{N}\|_{\mathrm{F}}/\|\mathbf{K}_{N}\|_{\mathrm{F}}$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
