clear
close all

% Add utils and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'utils');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'data_online');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'data_chapter7');
addpath(genpath(targetPath));

%%%%% CODE FOR FIBONACCI HAMILTNOIAN EXAMPLE (CASE STUDY 2) %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot the covers for various lambda

k = 10;
lamvec = 0:0.01:4;

Spec = cell(1,length(lamvec));
Spec2d = cell(1,length(lamvec));
Spec3d = cell(1,length(lamvec));

for jj = 1:length(lamvec)
    S1 = fib_siglam(k,lamvec(jj),'A');
    S2 = fib_siglam(k+1,lamvec(jj),'A');
    Spec{jj} = set_union(S1,S2,0);
    Spec2d{jj} = set_sum_lowmem_B(Spec{jj},Spec{jj},0);
    Spec3d{jj} = set_sum_lowmem_B(Spec2d{jj},Spec{jj},0);
end


follicle_plot(Spec, lamvec)
title('$S_{\lambda,k}, k=10$','interpreter','latex')
ylabel('$\lambda$','interpreter','latex')
ax=gca; ax.FontSize=18;
set(gcf, 'Position', [300 98 560*1.5 420]);

follicle_plot(Spec2d, lamvec)
title('$S_{\lambda,k}+S_{\lambda,k}, k=10$','interpreter','latex')
ylabel('$\lambda$','interpreter','latex')
ax=gca; ax.FontSize=18;
set(gcf, 'Position', [300 98 560*1.5 420]);

follicle_plot(Spec3d, lamvec)
title('$S_{\lambda,k}+S_{\lambda,k}+S_{\lambda,k}, k=10$','interpreter','latex')
ylabel('$\lambda$','interpreter','latex')
ax=gca; ax.FontSize=18;
set(gcf, 'Position', [300 98 560*1.5 420]);


%% Compute the L (bounds on error in Hausdorff metric)

kvec = 5:5:15;
lamvec = 0.1:0.1:10;
L = zeros(kvec(end),length(lamvec));

for k = kvec
    for jj = 1:length(lamvec)
        jj
        S1 = fib_siglam(k,lamvec(jj),'A');
        S2 = fib_siglam(k+1,lamvec(jj),'A');
        L(k,jj) = max(max(S1(:,2)-S1(:,1)),max(S2(:,2)-S2(:,1)));
    end
end

CC = gray(5); CC(end,:) = []; CC = flipud(CC);

figure
semilogy(lamvec,L(5,:),'linewidth',1,'color',CC(1,:))
hold on
semilogy(lamvec,L(10,:),'linewidth',1.5,'color',CC(2,:))
semilogy(lamvec,L(15,:),'linewidth',2,'color',CC(3,:))
xlabel('$\lambda$','interpreter','latex','fontsize',18)
title('$L_{\lambda,k}$','interpreter','latex','fontsize',18)
legend({'$k=5$','$k=10$','$k=15$'},'fontsize',16,'interpreter','latex','location','southwest')
ax=gca; ax.FontSize=18;
grid minor

%% Compute the box-counting dimension

clear
load('fib_1D.mat') % load spectral computations - these a precomputed spectra using the above method, book uses larger k for smaller lambda

%%
DELTA = 2.^(-10:-0.5:-100);
DELTA(DELTA>10^(-10))
SamS = 20;
dimB = zeros(length(lamvec),SamS) + NaN;

for jj=1:length(lamvec)
    jj
	s1 = Spec{20,jj};
    if lamvec(jj)>5
        DELTA = DELTA(DELTA<10^(-4));
    end
    dd = DELTA(DELTA>Er_bound(20,jj)); % adaptive delta depending on error bound
    for jjj=1:SamS
	    v = box_dim(double(s1),dd);
	    p=polyfit(log(dd),log(exp(v.*abs(log(dd)))),1);
        dimB(jj,jjj) = abs(p(1));
    end
end

dimB = mean(dimB,2);

%%

f=figure;
plot(lamvec,dimB,'linewidth',1,'color','k')
hold on
ylim([0,1])
ax=gca; ax.FontSize=16;
xlabel('$\lambda$','interpreter','latex')
grid minor
II=(lamvec>=4);
III=(lamvec>=8);
SL=2*lamvec+22;
SU=(lamvec-4+sqrt((lamvec-4).^2-12))/2;
plot(lamvec(III),log(1+sqrt(2))./log(SU(III)),'--k','linewidth',0.5)
plot(lamvec(II),log(1+sqrt(2))./log(SL(II)),'--k','linewidth',0.5)
f.Position=[360 173.6667 700 344];
grid minor

%% Compute number of connected components for cubic Fibonacci

k = 12;
lamvec = 1.5:0.02:3;
CC = zeros(1,length(lamvec));

for jj = 1:length(lamvec)
    S1 = fib_siglam(k,lamvec(jj),'A');
    S2 = fib_siglam(k+1,lamvec(jj),'A');
    Spec = set_union(S1,S2,0);
    Spec2d = set_sum_lowmem_B(Spec,Spec,0);
    Spec3d = set_sum_lowmem_B(Spec2d,Spec,0);
    CC(jj) = size(Spec3d,1);
end

figure
plot(lamvec,CC,'k.','markersize',12)
xlabel('$\lambda$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=16;


        