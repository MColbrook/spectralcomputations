clear
close all

% Add utils and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'utils');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'data_online');
addpath(genpath(targetPath));

%%%%%%%% CODE FOR EIGENVALUE SCHRODINGER OPERATOR EXAMPLE %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set parameters

phi = @(x) sin(x); % test function
u0 = @(x) x.^2./(1+x.^6)*sqrt(9/pi); % vector for scalar-valued spectral measure
V = @(x) -(3/2)./(1+x.^2); % potential
cut_off = 10^(-30); % cut-off parameter to make operator banded
SCALE = 1/10; % scale factor for basis functions
Nvec = unique(round((10.^(1:0.1:4))/2+1)); % vector of truncation sizes
Nvec = Nvec(Nvec<3500/2); % alterantive method than direct diagonalisation recommended for larger N

%% compute approximations of spectral measures

Int = zeros(length(Nvec),1);
for jjj=1:length(Nvec)

    N = Nvec(jjj) % truncation size is 2N+1

    % derivative matrix in terms of basis
    N=N+15; % trick to get actual truncations of matrices
    D1=spdiags(transpose([-N:N+2; (-(N+1):(N+1))+1/2]),[-1,0],2*N+3,2*N+3);
    D1=1i/2*(D1+D1')*SCALE;
    D2=D1*D1;
    D2=D2(2:end-1,2:end-1);
    N=N-15;
    DIFF=D2(16:end-15,16:end-15);

    % compute the Fourier series of V after mapping to unit circle
    v=chebfun(@(x) V(1i*(1-exp(1i*x))./(SCALE*(1+exp(1i*x)))) ,[-pi,pi],'trig');
    c=trigcoeffs(v,4*N+1);
    c(abs(c)<cut_off)=0;
    V1=sptoeplitz(c(2*N+1:4*N+1),c(2*N+1:(-1):1));

    A = -DIFF+V1; A = (A+A')/2;
    [VV,E]=eigs(A,size(A,1)); E = diag(E);

    f0=chebfun(@(x) u0(1i*(1-exp(1i*x))./(SCALE*(1+exp(1i*x))))./(1+exp(1i*x)) ,[-pi,pi],'trig');
    f0=trigcoeffs(f0,2*N+1)*(2*sqrt(pi/SCALE));
    f0=f0/norm(f0);
    c = abs(VV'*f0).^2;

    [E,I]=sort(E,'ascend');
    c = c(I);
    Int(jjj) = sum(phi(E(:)).*c(:)); % approximation of integral
end

%% Plot the results

figure
loglog(Nvec(1:end-1)*2+1,abs(Int(1:end-1)-Int(end)),'.-k','linewidth',1,'markersize',12)
xlabel('$n$','interpreter','latex','fontsize',18)
title('$\left|\int\sin(y)\mathrm{d}\mu_v(y)-\int\sin(y)\mathrm{d}\mu_{n;v}(y)\right|$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;


