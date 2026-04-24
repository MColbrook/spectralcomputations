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


%%%%%%%%%%%%%%% CODE FOR THE ALMOST MATHIEU EXAMPLE %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% EXPERIMENTS IN CHAPTER INTRO %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Draw Hofstadter butterfly

qmax = 101;
alf = [];
for q=2:qmax 
    p = (1:q-1)';
    p = p(gcd(p,q)==1);
    alf = [alf;[p/q p q*ones(length(p),1)]];
end 
[~,indx] = sort(alf(:,1));
alf = alf(indx,:);

f=figure;
clf

pf = parfor_progress(size(alf,1));
pfcleanup = onCleanup(@() delete(pf));

for j=1:size(alf,1)
    alpha = alf(j,1);
    p = alf(j,2);
    q = alf(j,3);
    r = gcd(p,q);
    p = round(p/r);
    q = round(q/r);
    S = AM_spec(p,q,1,0.003);

    bw = diff(S');                % test bandwidths for plot
    if min(bw)<=0
        fprintf('warning: alpha = %10.7f (p=%3d; q=%3d) has %d bands have non-positive width:\n', alpha, p, q, sum(bw<=0));
        indx = find(bw<=0);
        for k=1:length(indx)
            fprintf('   [%20.15f,%20.15f] has width %10.5e\n',S(indx(k),1),S(indx(k),2),bw(indx(k)))
        end
    end
    S = reshape([S NaN*ones(size(S,1),1)]',3*size(S,1),1)+1i*alpha*ones(3*size(S,1),1);
    plot(real(S),imag(S),'-','linewidth',0.12*2,'color','k')
    hold on
    parfor_progress(pf);
end

title('$\mathrm{Sp}_{+}(\alpha,1)$','interpreter','latex')
ylabel('$\alpha$','interpreter','latex')
ax=gca; ax.FontSize=18;

f.Position=[246.3333   79.6667  890.6667  522.6667];

%%%%%%%%%%%%%%%%%% EXPERIMENTS IN Section 7.2.5 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Draw lambda dependent spectrum

q=300;
sigma = cell(q+1,1);
for p=0:q
    sigma{p+1} = AM_spec(377,610,p/q,0);
end

figure
pm_follicle_AM(sigma,q)
ylim([0,1])
title('$\mathrm{Sp}_{+}(377/610,\lambda)$','interpreter','latex')
ylabel('$\lambda$','interpreter','latex')
ax=gca; ax.FontSize=18;


%% Compute quantities that are continuous under S1

lamvec=0:0.05:1;
pvec = [3,5,8,13,21,34,55,89];%,144,233,377,610,987,1597,2584,4181,6765,10946];
for k=1:11
    pvec=[pvec,sum(pvec(end-1:end))];
end
qvec = pvec(2:end);
pvec = pvec(5:5:end);
qvec = qvec(5:5:end);

er = zeros(length(qvec),length(lamvec));
sig = cell(length(qvec),length(lamvec));

for jj=1:length(qvec)
    pvec(jj)
    pf = parfor_progress(length(lamvec));
    pfcleanup = onCleanup(@() delete(pf));
    for ii=1:length(lamvec)
        er(jj,ii) = 13.2*sqrt(lamvec(ii))*sqrt(abs(pvec(jj)/qvec(jj)-(-1+sqrt(5))/2)); 
        sig{jj,ii} = AM_spec(pvec(jj),qvec(jj),lamvec(ii),er(jj,ii));
        vv = sig{jj,ii};
        parfor_progress(pf);
    end
end

%% Lebesgue measure

CC = gray(6); CC(end,:) = []; CC = flipud(CC); % for greyscale plot

leb = 0*er;
for jj=1:3
    for ii=1:size(er,2)
        vv=sig{jj,ii};
        leb(jj,ii)=sum(vv(:,2)-vv(:,1));
    end
end

figure
plot(lamvec,leb(1,:),'.-','linewidth',1,'markersize',12,'color',CC(1,:))
hold on
semilogy(lamvec,leb(2,:),'*-','linewidth',1,'markersize',6,'color',CC(2,:))
semilogy(lamvec,leb(3,:),'+-','linewidth',1,'markersize',5,'color',CC(3,:))
grid minor
xlim([0,1]); ylim([0,4.5])

title('Lebesgue Measure','interpreter','latex')
xlabel('$\lambda$','interpreter','latex')
ax=gca; ax.FontSize=24;
hleg = legend({sprintf('$%d/%d$',pvec(1),qvec(1)),sprintf('$%d/%d$',pvec(2),qvec(2)),...
    sprintf('$%d/%d$',pvec(3),qvec(3))},...
    'interpreter','latex', 'fontsize',14,'location','best');
title(hleg,'$p/q$','interpreter','latex', 'fontsize',14)


%% Number of components

CC = gray(6); CC(end,:) = []; CC = flipud(CC);

comp = 0*er;
for jj=1:size(er,1)
    for ii=1:size(er,2)
        comp(jj,ii)=size(sig{jj,ii},1);
    end
end

figure
semilogy(lamvec,comp(1,:),'.-','linewidth',1,'markersize',12,'color',CC(1,:))
hold on
semilogy(lamvec,comp(2,:),'*-','linewidth',1,'markersize',6,'color',CC(2,:))
semilogy(lamvec,comp(3,:),'+-','linewidth',1,'markersize',5,'color',CC(3,:))
grid minor
xlim([0,1])
title('\# Connected Components','interpreter','latex')
xlabel('$\lambda$','interpreter','latex')
ax=gca; ax.FontSize=24;
hleg = legend({sprintf('$%d/%d$',pvec(1),qvec(1)),sprintf('$%d/%d$',pvec(2),qvec(2)),...
    sprintf('$%d/%d$',pvec(3),qvec(3))},...
    'interpreter','latex', 'fontsize',14,'location','best');
title(hleg,'$p/q$','interpreter','latex', 'fontsize',14)


%% Capacity

for jj=1:length(qvec)
    pvec(jj)
    pf = parfor_progress(length(lamvec));
    pfcleanup = onCleanup(@() delete(pf));
    for ii=1:length(lamvec)
        cap(jj,ii) = spec_capacity(sig{jj,ii});
        parfor_progress(pf);
    end
end

CC = gray(6); CC(end,:) = []; CC = flipud(CC);
figure
plot(lamvec,cap(1,:),'.-','linewidth',1,'markersize',12,'color',CC(1,:))
hold on
semilogy(lamvec,cap(2,:),'*-','linewidth',1,'markersize',6,'color',CC(2,:))
semilogy(lamvec,cap(3,:),'+-','linewidth',1,'markersize',5,'color',CC(3,:))
grid minor
axis([0,1,0.98,1.4])
title('Capacity','interpreter','latex')
xlabel('$\lambda$','interpreter','latex')
ax=gca; ax.FontSize=24;
hleg = legend({sprintf('$%d/%d$',pvec(1),qvec(1)),sprintf('$%d/%d$',pvec(2),qvec(2)),...
    sprintf('$%d/%d$',pvec(3),qvec(3))},...
    'interpreter','latex', 'fontsize',14,'location','best');
title(hleg,'$p/q$','interpreter','latex', 'fontsize',14)


%% Box-counting dimension

clear
% load('AM_critical_golden.mat') % high-res computed spectra
load('AM_critical_Cahen.mat') % high-res computed spectra

S = sig{end};
er = er(end);
M = ceil(log2(1/10^(-9)*2));
delta = 2.^(-(5:0.05:M));
delta = delta(delta<10^(-2));
SamS = 10; % number of samples for random end point perturbation (results may vary for large delta depending on seed)
dim = zeros(SamS,length(delta));

for j=1:SamS
    dim(j,:) = box_dim(S,delta,0);
end

NN = exp(log(1./delta).*mean(dim,1));
NN = NN/NN(1);

c = log(NN)./log(delta(1)./delta);

%% Plot the box-counting dimension results

figure
semilogx(delta(delta>er),c(delta>er),'k','linewidth',2)
hold on
plot(delta(delta>er),0.5+0.*delta(delta>er),'--k','linewidth',1)
xlabel('$\delta$','interpreter','latex')
ax=gca; ax.FontSize=18;
title('$\log(N_{\delta}(\mathrm{Sp}_{+}(\alpha,1))/N_{\delta_0}(\mathrm{Sp}_{+}(\alpha,1)))/\log(\delta_0/\delta)$','interpreter','latex','fontsize',16)
xlim([min(delta(delta > er)),0.01])
ylim([0.48,0.7])
grid minor

figure
loglog(delta,exp(log(1./delta).*mean(dim,1)),'k.-','linewidth',1,'markersize',15)
hold on
XX = delta(2:end-3);
loglog(XX,6.5*XX.^(-1/2),':k','linewidth',2)
loglog(XX,4.5*XX.^(-2/3),':k','linewidth',2)
ax=gca; ax.FontSize=18;
xlabel('$\delta$','interpreter','latex')
title('$N_{\delta}(\mathrm{Sp}_{+}(\alpha,1))$','interpreter','latex')
xlim([10^(-10),.01])
ylim(10.^([2,7]))
grid minor
text(10^(-7),10000/1.1,'$\mathcal{O}(1/\delta^{1/2})$','interpreter','latex','fontsize',18,'Rotation',-32)
text(10^(-7),1000000/2,'$\mathcal{O}(1/\delta^{2/3})$','interpreter','latex','fontsize',18,'Rotation',-37)

%% Hausdorff dimension

clear
load('AM_critical_golden.mat')

dimH = Haus_dim(sig{end-1},5:28,5:5:15,10^(-10),3);

CC = gray(3); CC(end,:) = CC(end,:)*0.9;
CC = flipud(CC);

figure
p1 = plot(5:28,dimH(:,1),'-^','markersize',7,'color',CC(1,:),'linewidth',1);
p1.MarkerEdgeColor = [0 0 0];
p1.MarkerFaceColor = CC(1,:);
hold on
p2 = plot(5:28,dimH(:,2),'-s','markersize',7,'color',CC(2,:),'linewidth',1);
p3 = plot(5:28,dimH(:,3),'-d','markersize',7,'color',CC(3,:),'linewidth',1);
p2.MarkerEdgeColor = [0 0 0];
p2.MarkerFaceColor = CC(2,:);
p3.MarkerEdgeColor = [0 0 0];
p3.MarkerFaceColor = CC(3,:);

ax=gca; ax.FontSize=18;
xlabel('$n_1$','interpreter','latex')
xlim([0,30])
ylim([0.3,0.7])
title('$\alpha=(\sqrt{5}-1)/2$','interpreter','latex')
hleg = legend({'$n_2=5$','$n_2=10$','$n_2=15$'},...
    'interpreter','latex', 'fontsize',14,'location','best');
grid minor




