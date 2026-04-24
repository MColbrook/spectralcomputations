clear
close all

% Add Koopman algorithms and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'Koopman_algorithms');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'Koopman_datasets');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%%%%%%%% CODE FOR FIGURE 11.28 %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1) % set random seed for reproducibility

%% Choose observable
g = @(x) abs(x-1/3)+(x>0.78)+sin(20*x); gg = chebfun(@(x)  g(x),[0,1],'splitting','on');
f = @(x) g(x)/sqrt(sum(gg*conj(gg))); ff = chebfun(@(x)  f(x),[0,1],'splitting','on');  % function we compute measure wrt
n=15; % trunction size to compute the inner products
K = Koopman_Haar(n); % koopman matrix wrt Haar system

%% Compute wavelet coefficients
xpts=1/(2^(n+1)):2^(-n):1;
[coeffs,L]=wavedec(f(xpts),round(log2(length(xpts))),'db1');
coeffs=coeffs/sqrt(length(xpts));
at=coeffs(1);
coeffs(1)=0; % subtract off the atomic part, will add on later
coeffs=coeffs(:);

%% Compute the fourier coefficients
MU=zeros(2*(2^n)+1,1);  N=500;
MU(N+1)=(coeffs'*coeffs)/(2*pi);

Fc=coeffs;
for j=1:N
    Fc=K*Fc;
    MU(N-j+1)=(coeffs'*Fc)/(2*pi);
    MU(N+1+j)=(MU(N-j+1))';
end

MU=MU+abs(at)^2/(2*pi); % add atomic part back on

%% Spectral measures plot
N1=100; N2=500;

MU01=chebfun(full(MU((N+1-N1):(N+1+N1))),[-pi pi],'trig','coeffs'); % approximation with no filter
MU1=MomentMeas(MU((N+1-N1):(N+1+N1)),'filt','vand');  % approximation with filter
MU02=chebfun(full(MU((N+1-N2):(N+1+N2))),[-pi pi],'trig','coeffs'); % approximation with no filter
MU2=MomentMeas(MU((N+1-N2):(N+1+N2)),'filt','vand');  % approximation with filter

%% N=100 plot
figure
semilogy(max(real(MU01),0),'color',[1,1,1]*0.5,'linewidth',0.5)
hold on
semilogy(max(real(MU1),10^(-16)),'color','k','linewidth',1);
xlabel('$\theta$','interpreter','latex','fontsize',18)
title('$N=100$','interpreter','latex','fontsize',18)
ylim([10^(-2),100]); ax = gca; ax.FontSize = 18;
legend({'no filter','with filter'},'interpreter','latex','fontsize',16,'location','northwest')

%% N=500 plot
figure
semilogy(max(real(MU02),0),'color',[1,1,1]*0.5,'linewidth',0.5)
hold on
semilogy(max(real(MU2),10^(-16)),'color','k','linewidth',1);
xlabel('$\theta$','interpreter','latex','fontsize',18)
title('$N=500$','interpreter','latex','fontsize',18)
ylim([10^(-2),100]); ax = gca; ax.FontSize = 18;
legend({'no filter','with filter'},'interpreter','latex','fontsize',16,'location','northeast')

a2 = axes();    a2.Position = [0.18 0.65 0.3 0.2]; % xlocation, ylocation, xsize, ysize
semilogy(max(real(MU02),0),'color',[1,1,1]*0.5,'linewidth',0.5)
hold on
semilogy(max(real(MU2),10^(-16)),'color','k','linewidth',1);
axis([-1,-0.8,10^(-5),1])



function K = Koopman_Haar(n)
K=sparse(2^n,2^n);
% set up the indexing
nId=0;    kId=0;
for j=1:n
    nId=[nId(:);    zeros(2^(j-1),1)+j];    kId=[kId(:);    transpose(0:(2^(j-1)-1))];
end
K(1,1)=1;
for j=2:2^(n-1)
    I1 = find((nId==(nId(j)+1)).*(kId==kId(j)));
    if ~isempty(I1)
        K(I1,j)=1/sqrt(2);
    end
    I1 = find((nId==(nId(j)+1)).*(kId==(2^(nId(j))-(kId(j)+1))));
    if ~isempty(I1)
        K(I1,j)=-1/sqrt(2);
    end
end
end
    