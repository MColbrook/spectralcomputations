clear
close all

% Add utils and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'utils');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'data_online');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%% CODE FOR Example in Section 6.4.1 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set parameters (change these to generate the different plots in the book)

c1 = 0.5;
c2vec = [1/4,3/4];
Phi = (sqrt(5)-1)/2;

n2vec = [5,10,20];
n1vec = round(10.^(1:0.1:5));

N = 2*max(n1vec)+100;

c = zeros(2*(N+2)+1,1);
c(N+3) = 1;

MUac = zeros(length(n1vec),length(n2vec),length(c2vec));
MUp = zeros(length(n1vec),length(n2vec),length(c2vec));

%% Perform computation

for ll = 1:length(c2vec)
    c2 = c2vec(ll);

    % Set up matrix

    al = @(n) (mod(n,2)==0).*c2.*sin(2*pi*(n/2*Phi)) + (mod(n,2)==1)*sqrt(1-c1^2);
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

    % Compute spectral type

    [mup,muac] = mu_type(U,c,n1vec,n2vec,n2vec);
    
    MUp(:,:,ll) = mup;
    MUac(:,:,ll) = muac;
end

%% Plot the results

cc=0.25;

figure
loglog(n1vec,MUac(:,n2vec==5,c2vec==cc),'-k','linewidth',2)
hold on
loglog(n1vec,MUac(:,n2vec==10,c2vec==cc),'--k','linewidth',2)
loglog(n1vec,MUac(:,n2vec==20,c2vec==cc),':k','linewidth',2)

legend({'$n_2=5$','$n_2=10$','$n_2=20$'},'fontsize',16,'interpreter','latex','location','southwest')
xlim([10,10^5])

xlabel('$n_1$','interpreter','latex','fontsize',18)
title('$c_2=1/4$','interpreter','latex','fontsize',18)

ax=gca; ax.FontSize=18;



figure
loglog(n1vec,MUp(:,n2vec==5,c2vec==cc),'-k','linewidth',2)
hold on
loglog(n1vec,MUp(:,n2vec==10,c2vec==cc),'--k','linewidth',2)
loglog(n1vec,MUp(:,n2vec==20,c2vec==cc),':k','linewidth',2)

legend({'$n_2=5$','$n_2=10$','$n_2=20$'},'fontsize',16,'interpreter','latex','location','southwest')
xlim([10,10^5])

xlabel('$n_1$','interpreter','latex','fontsize',18)
title('$c_2=1/4$','interpreter','latex','fontsize',18)

ax=gca; ax.FontSize=18;




