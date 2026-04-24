clear
close all

% Add Koopman algorithms and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'Koopman_algorithms');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'Koopman_datasets');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%%%%%%%% CODE FOR FIGURE 11.27 %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1) % set random seed for reproducibility

%% Set parameters for numerical integration of ODE
options = odeset('RelTol',1e-14,'AbsTol',1e-14);
SIGMA = 10;   BETA = 8/3;   RHO = 28;
ODEFUN = @(t,y) [SIGMA*(y(2)-y(1));y(1).*(RHO-y(3))-y(2);y(1).*y(2)-BETA*y(3)];

%% Set parameters for the example
M = 10^6;    % number of snapshots      
dt = 0.05;    % time step for trajectory sampling

%% Produce the data
Y0 = (rand(3,1)-0.5)*4; % initial point off the Lorenz attractor
[~,Y0] = ode45(ODEFUN,[0.000001 1, 100],Y0,options); Y0 = Y0(end,:)'; % new initial point on the Lorenz attractor
[~,DATA] = ode45(ODEFUN,[0.000001 (1:(M-1))*dt],Y0,options);

%% Compute moments and measures
g = @(x,y,z) tanh((x.*y-3*z)/5);
X = g(DATA(:,1),DATA(:,2),DATA(:,3));
c1 = sqrt(mean(X.*X));
X = X/c1;
c2 = mean(X);
X = X - c2;

MU = ErgodicMoments(X(1:10^4),1000);
mu1 = MomentMeas(MU,'filt','vand');
MU = ErgodicMoments(X(1:10^6),1000);
mu2 = MomentMeas(MU,'filt','vand');


figure
semilogy(mu1,'color',[1,1,1]*0.7,'linewidth',2)
hold on
semilogy(mu2,'color','k','linewidth',1);
xlim([0,pi])
ylim([0.001,100])
title('$\xi_{N;g}^{\sigma}(\theta),N=10^3$','interpreter','latex','fontsize',18)
xlabel('$\theta$','interpreter','latex','fontsize',18)
legend({'$M=10^4$','$M=10^6$'},'interpreter','latex','fontsize',16,'location','northeast')
ax = gca; ax.FontSize = 18;


%% Convergence of Fourier coefficient wrt M

Mvec = unique(round(10.^([2.01:0.01:5,5.3])));
mu = 0*Mvec;

for jj = 1: length(Mvec)
    MU = ErgodicMoments(X(1:Mvec(jj)),100);
    mu(jj) = MU(end);
end


figure
loglog(Mvec,abs(mu-mu(end)),'k','linewidth',1)
hold on
loglog(Mvec(2:end-5),Mvec(2:end-5).^(-1/2),'--k','linewidth',1)
xlim([100,10^5])
ylim([10^(-5)*2,0.1])
title('Error of $\widehat{\xi_{g}}(100)$','interpreter','latex','fontsize',18)
text(10000/3,0.03,'$\mathcal{O}(M^{-1/2})$','interpreter','latex','fontsize',18,'Rotation',-17)
xlabel('$M$','interpreter','latex','fontsize',18)
ax = gca; ax.FontSize = 17;
