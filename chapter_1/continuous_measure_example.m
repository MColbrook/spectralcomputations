clear
close all

% Add utils and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'utils');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'data_online');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%% CODE FOR EXAMPLE IN SECTION 1.1.3 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set parameters
u0 = @(x) x.^2./(1+x.^6)*sqrt(9/pi); % function we compute measure wrt
V = @(x) -(3/2)./(1+x.^2); % potential
cut_off = 10^(-30); % cut-off parameter to make operator banded
N = 50; % size of discretization is 2N+1
SCALE=1/10; % scale factor

%% Derivative matrix in terms of basis
N=N+15; % trick to get actual truncations of matrices
D1=spdiags(transpose([-N:N+2; (-(N+1):(N+1))+1/2]),[-1,0],2*N+3,2*N+3);
D1=1i/2*(D1+D1')*SCALE;
D2=D1*D1;
D2=D2(2:end-1,2:end-1);
N=N-15;
DIFF=D2(16:end-15,16:end-15);

%% Compute the Fourier series of V after mapping to unit circle
v=chebfun(@(x) V(1i*(1-exp(1i*x))./(SCALE*(1+exp(1i*x)))) ,[-pi,pi],'trig');
c=trigcoeffs(v,4*N+1);
c(abs(c)<cut_off)=0;
V1=sptoeplitz(c(2*N+1:4*N+1),c(2*N+1:(-1):1));

A=-DIFF+V1;
A=(A+A')/2;

[V,E]=eig(full(A),'vector');

f0=chebfun(@(x) u0(1i*(1-exp(1i*x))./(SCALE*(1+exp(1i*x))))./(1+exp(1i*x)) ,[-pi,pi],'trig');
f0=trigcoeffs(f0,2*N+1)*(2*sqrt(pi/SCALE));
f0=f0/norm(f0);

c = abs(V'*f0).^2; % amplitudes for measure

[E,I]=sort(E,'ascend');
c = c(I);
V = V(:,I);
E3=E;

E = [E(:),c(:)*1i+E(:),0*E(:)+NaN]; % trick for MATLAB for bar plot
E=conj(E'); E=E(:);

%% Plot the measure

ff=figure;
semilogy(real(E),imag(E)+10^(-10),'k','linewidth',0.8);
xlim([-1,2])
ylim([10^(-6),1])
load('cts_spec.mat')
hold on
mu1(1:2)=mu1(3);
mu1(1:5)=smooth(mu1(1:5));
plot(X,mu1,'k','linewidth',2)
title(sprintf('Spectral measure $\\mu_{V_N}$ ($N=%d$)',2*N+1),'fontsize',12,'interpreter','latex');
xlabel('Spectral parameter $\lambda$','interpreter','latex','fontsize',12)
ax=gca; ax.FontSize=12;
ff.Position=[525.6667  343.0000  560.0000  284.6667];

%% Plot the cdf

TH = E(1:3:end);
[~,Ib] = sort(TH(:),'ascend');
THp=TH(Ib); THp=[THp(:)-10^(-14),THp(:)]'; THp=THp(:);
cdf=0*THp+10^(-14);
cc=10^(-14);
for j=1:length(TH)
    cdf(2*j-1)=cc;    cc=cc+c(Ib(j));    cdf(2*j)=cc;
end

THp = [-1;THp(:);1000000]; cdf = [10^(-14);cdf(:);sum(c)]; % for visualisation


ff=figure;
plot(THp,cdf,'k','linewidth',0.8)
xlim([-1,2])
ylim([0.61,0.63])
title(sprintf('Cdf ($N=%d$)',2*N+1),'interpreter','latex','fontsize',14)
xlabel('Spectral parameter $\lambda$','interpreter','latex','fontsize',14)
ax=gca; ax.FontSize=14;
ff.Position=[525.6667  343.0000  560.0000  284.6667];

%% Plot zoom in for all three N in book (adjust above code to produce data)

clear
load('cdf101.mat')
load('cdf1001.mat')
load('cdf10001.mat')
ff=figure;
plot(THp1,cdf1,'color',[1,1,1]*0.8,'linewidth',3.2)
hold on
plot(THp2,cdf2,'color',[1,1,1]*0.4,'linewidth',1.8)
plot(THp3,cdf3,'color',[1,1,1]*0,'linewidth',0.8)
xlim([-1,2])
ylim([0.615,0.63])
xlabel('Spectral parameter $\lambda$','interpreter','latex','fontsize',14)
ax=gca; ax.FontSize=14;
ff.Position=[525.6667  343.0000  560.0000  284.6667];
legend({'$N=101$','$N=1001$','$N=10001$'},'fontsize',14,'interpreter','latex','location','northwest')


a2 = axes();
a2.Position = [0.67 0.3 0.2 0.3]; % xlocation, ylocation, xsize, ysize
plot(THp1,cdf1,'color',[1,1,1]*0.8,'linewidth',3.2)
hold on
plot(THp2,cdf2,'color',[1,1,1]*0.4,'linewidth',1.6)
plot(THp3,cdf3,'color',[1,1,1]*0,'linewidth',0.8)
axis([0.1,0.3,0.62,0.6225])

