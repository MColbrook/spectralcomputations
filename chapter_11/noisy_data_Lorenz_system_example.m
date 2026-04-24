clear
close all

% Add Koopman algorithms and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'Koopman_algorithms');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'Koopman_datasets');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%%%%%%% CODE FOR FIGURE 11.26 %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1) % set random seed for reproducibility

%% Set parameters for numerical integration of ODE
options = odeset('RelTol',1e-14,'AbsTol',1e-14);
SIGMA = 10;   BETA = 8/3;   RHO = 28;
ODEFUN = @(t,y) [SIGMA*(y(2)-y(1));y(1).*(RHO-y(3))-y(2);y(1).*y(2)-BETA*y(3)];

%% Set parameters for the example
N = 200;    % number of delays (this is an upper bound, below we reduce at various stages)
M = 10^5;    % number of snapshots      
dt = 0.05;    % time step for trajectory sampling

%% Produce the data
M2 = M; % for visualisation of functions we use more data points
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

%% Produce the data for the error curves
N = 50;
g = @(x,y,z) tanh((x.*y-3*z)/5);    % observable
nvec = 0:0.02:1; % noise levels
sam = 100; % number of samples (reduce for speed if wanted)

for M = [1000,10000,100000]
    R = zeros(2,length(nvec),sam);
    
    PX = g(PX1(1:M,1:N),PX2(1:M,1:N),PX3(1:M,1:N));
    PY = g(PY1(1:M,1:N),PY2(1:M,1:N),PY3(1:M,1:N));
    
    cc = sqrt(norm(PX(:,1))^2/M);
    PX = PX/cc;
    PY = PY/cc;
    
    for jj = 1:length(nvec)
        jj
        for ii = 1:sam
            Nr = nvec(jj)*mean(std(PX))*randn(size(PX,1),N);
            Nr2 = nvec(jj)*mean(std(PX))*randn(size(PX,1),N);
            
            PX0 = PX + Nr; % add noise to measurements
            PY0 = PY + Nr2; % add noise to measurements
            
            K = PX0\PY0;
            [V,LAM] = eig(K,'vector');
            R(1,jj,ii) = mean(vecnorm(PY*V-PX*V*diag(LAM))./vecnorm(PX*V));
            
            [~,V,mpD] = mpEDMDqr(PX0,PY0,1/M);
            LAM = diag(mpD);
            R(2,jj,ii) = mean(vecnorm(PY*V-PX*V*diag(LAM))./vecnorm(PX*V));
        end
    end

    R = squeeze(mean(R,3));
    save(sprintf('noisy_out%d.mat',round(log10(M)-2)),'nvec','R')
end


%% Plot the results
clear
close all
nvec = 0:0.02:1;

figure
load('noisy_out1.mat')
plot(nvec,R(1,:),'k','linewidth',2)
hold on
load('noisy_out2.mat')
plot(nvec,R(1,:),'--k','linewidth',1)
load('noisy_out3.mat')
plot(nvec,R(1,:),':k','linewidth',1)
ylim([0,0.8])
grid minor
legend({'$M=10^3$','$M=10^4$','$M=10^5$'},'fontsize',16,'interpreter','latex','location','northwest')
xlabel('Noise level','interpreter','latex','fontsize',18)
title('EDMD Eigenpair Residuals','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

figure
load('noisy_out1.mat')
plot(nvec,R(2,:),'k','linewidth',2)
hold on
load('noisy_out2.mat')
plot(nvec,R(2,:),'--k','linewidth',1)
load('noisy_out3.mat')
plot(nvec,R(2,:),':k','linewidth',1)
ylim([0,0.8])
grid minor
legend({'$M=10^3$','$M=10^4$','$M=10^5$'},'fontsize',16,'interpreter','latex','location','northwest')
xlabel('Noise level','interpreter','latex','fontsize',18)
title('mpEDMD  Eigenpair Residuals','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;




