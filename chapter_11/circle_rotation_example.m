clear
close all

% Add Koopman algorithms and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'Koopman_algorithms');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'Koopman_datasets');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%%%%%%%% CODE FOR FIGURE 11.30 %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1D example for the point spectrum

TRUNC=1000;
a = 2*pi*0.7;
% a = 2*pi/sqrt(2);

% compute Fourier expansion of g
g = chebfun(@(x)  cos(11*x)/(sqrt(1.001+sin(x))) ,[-pi,pi],'trig');
coeffs = trigcoeffs(g,2*TRUNC+1);
coeffs(abs(coeffs)<10^(-14))=0;

% compute filtered measures
N=100;
MU = zeros(2*N+1,1)';
Fc=coeffs;
MU(N+1)=(coeffs'*coeffs);
for j=1:N
    K = exp((-TRUNC:TRUNC)'*1i*a*j);
    Fc=K.*coeffs;
    MU(N+1+j)=(Fc'*coeffs);
    MU(N-j+1)=(MU(N+1+j))';
end
mu1=MomentMeas(MU,'filt','fejer');  % approximation with filter
mu1=real(mu1)/sum(real(mu1));
K0=MomentMeas(MU*0+1,'filt','fejer');   K0=K0(0);
mu1 = mu1/real(K0);

N=10000;
MU = zeros(2*N+1,1)';
Fc=coeffs;
MU(N+1)=(coeffs'*coeffs);
for j=1:N
    K = exp((-TRUNC:TRUNC)'*1i*a*j);
    Fc=K.*coeffs;
    MU(N+1+j)=(Fc'*coeffs);
    MU(N-j+1)=(MU(N+1+j))';
end
mu2=MomentMeas(MU,'filt','fejer');  % approximation with filter
mu2=real(mu2)/sum(real(mu2));
K0=MomentMeas(MU*0+1,'filt','fejer');   K0=K0(0);
mu2 = mu2/real(K0);


% plot the results

figure
subplot(2,1,1)
semilogy(mu1,'k','linewidth',0.5)
ax = gca; ax.FontSize = 14;
ylim([10^(-10),1])
title('$N=100$','interpreter','latex','fontsize',14)

subplot(2,1,2)
semilogy(mu2,'k','linewidth',0.5)
ax = gca; ax.FontSize = 14;
ylim([10^(-10),1])
title('$N=10000$','interpreter','latex','fontsize',14)
xlabel('$\theta$','interpreter','latex','fontsize',18)

