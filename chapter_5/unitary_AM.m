clear
close all

% Add utils and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'utils');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'data_online');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%% CODE FOR Example in Section 5.5.2 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set model parameters

lam1 = 1/sqrt(2);
lam2 = 1/sqrt(2);
Theta = 0;
Phi_vec = 0:0.01:0.5;

%% Set algorithmic parameters - increase order, N and spacing of X or higher resolution

X = -pi:0.02:0;
L = length(X);
order = 2;
epsilon = 0.02;
[poles,res]=rational_kernel(order,'equi');
N = 10000; % can also diagonalise polar decomposition for improved efficiency
c = zeros(2*(N+2)+1,1);
c(N+3) = 1;
mu = zeros(length(Phi_vec),L);

%% Compute butterfly

pf = parfor_progress(L*length(Phi_vec));
pfcleanup = onCleanup(@() delete(pf));

for jj=1:length(Phi_vec)

    Phi = Phi_vec(jj);

    % set up matrix

    al = @(n) (mod(n,2)==0).*lam2.*sin(2*pi*(n/2*Phi+Theta))    +    (mod(n,2)==1)*sqrt(1-lam1^2);
    rh = @(n) sqrt(1-abs(al(n)).^2);
    nn = (-(N+2):N+2)';

    % no conjugates needed for this example

    D0 = -al(nn).*al(nn-1);

    D11 = al(nn+1).*rh(nn); D11(N+2:2:end) = 0; D11(N:(-2):1) = 0;
    D12 = -al(nn-1).*rh(nn); D12(N+1:2:end) = 0; D12(N-1:(-2):1) = 0;
    D1 = D11+D12;
    D1 = [0;D1(1:end-1)];

    D2 = rh(nn).*rh(nn+1); D2(N+2:2:end) = 0; D2(N:(-2):1) = 0;
    D2 = [0;0;D2(1:end-2)];

    D11m = al(nn).*rh(nn-1); D11m(N+2:2:end) = 0; D11m(N:(-2):1) = 0;
    D12m = -al(nn-2).*rh(nn-1); D12m(N+1:2:end) = 0; D12m(N-1:(-2):1) = 0;
    D1m = D11m+D12m;
    D1m = [D1m(2:end);0];

    D2m = rh(nn-1).*rh(nn-2); D2m(N+1:2:end) = 0; D2m(N-1:(-2):1) = 0;
    D2m = [D2m(3:end);0;0];

    U = spdiags([D2m,D1m,D0,D1,D2],[-2,-1,0,1,2],2*(N+2)+1,2*(N+2)+1);

    % compute spectral measure

    for j=1:L
        for mm=1:order
            b = (U+exp(1i*X(j)-1i*epsilon*poles(mm))*speye(size(U)))'*c;
            a = res(mm)*((U-exp(1i*X(j)-1i*epsilon*poles(mm))*speye(size(U)))\c);
            mu(jj,j) = mu(jj,j)-real(b'*a)/pi/2;
        end
        parfor_progress(pf);
    end
end

%% Plot result

II = find(Phi_vec>0.01);
f=figure;
imagesc([X(1),-X(1)],[Phi_vec(II(1)),1-Phi_vec(II(1))],log10(max([mu(II,:),fliplr(mu(II,1:end-1));flipud([mu(II(1:end-1),:),fliplr(mu(II(1:end-1),1:end-1))])],10^(-16))));
clim([-6,2])
colormap turbo
set(gca,'ydir','normal');
axis tight
xlim([0,pi])
ylim([0,1])
xlabel('$\theta$','interpreter','latex','fontsize',18)
ylabel('$\Phi$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
f.Position=[360.0000   97.6667  550  420.0000];


