clear
close all
rng(0)

% Add AM_spec in folder from chapter 7 (computes spectra of AM operator)
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'chapter_7', 'routines');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%% CODE FOR THE ALMOST MATHIEU EXAMPLE %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set parameters

n1 = 1000;
n2vec = 0:2; % n2 parameters for the plot

lambda = 1; % coupling parameter
p = 2; q = 3; % alpha = p/q

%% Compute essential spectrum and distance function (Figure 8.3)

SPess = AM_spec(p,q,lambda,0); % spectrum without perturbation (essential spectrum)
xvec = -2:0.001:-1.5;

nn = -(n1+1):(n1+1);
pot = -rand(length(nn),1)*0+1./(abs(nn')+1); % random potential (generates discrete spectra)
H = spdiags([ones(2*n1+3,1),ones(2*n1+3,1),2*lambda*cos(2*pi*p/q*nn')+pot],[-1,1,0],2*n1+3,2*n1+3);
I = speye(size(H));
[V,D] = eig(full(H(2:end-1,2:end-1)));

nn = nn(2:end-1);
H = H(:,2:end-1);
I = I(:,2:end-1);
RES = sqrt(dot(H*V-I*D*V,H*V-I*V*D)); % use residual to avoid spurious eigenvalues

tau = zeros(length(n2vec),length(xvec));
pf = parfor_progress(length(xvec)*length(n2vec));
pfcleanup = onCleanup(@() delete(pf));
for ii=1:length(n2vec)
    for jj=1:length(xvec)
        Indx = find(abs(nn)>n2vec(ii)-1);
        tau(ii,jj) = svds(H(:,Indx)-xvec(jj)*I(:,Indx),1,'smallest');
        parfor_progress(pf);
    end
end

%% Plot the results

ff=figure;
E=diag(D);

subplot(3,1,1)
plot(xvec,tau(1,:),'k','linewidth',1)
hold on
plot(SPess(1,:)',0*SPess(1,:)','k','linewidth',4)
plot(E((RES'<0.0001)&(E>SPess(1,2)+0.001)),E((RES'<0.0001)&(E>SPess(1,2)+0.001))*0,'ok','linewidth',1,'markersize',8)
axis([-2.1,-1.6,0,0.15])
ax=gca; ax.FontSize=16;

text(-2.05,0.1,'$n_2=0$','interpreter','latex','fontsize',18,'Rotation',0)

subplot(3,1,2)
plot(xvec,tau(2,:),'k','linewidth',1)
hold on
plot(SPess',0*SPess','k','linewidth',4)
plot(E((RES'<0.0001)&(E>SPess(1,2)+0.001)),E((RES'<0.0001)&(E>SPess(1,2)+0.001))*0,'ok','linewidth',1,'markersize',8)
axis([-2.1,-1.6,0,0.15])
ax=gca; ax.FontSize=16;
text(-2.05,0.1,'$n_2=1$','interpreter','latex','fontsize',18,'Rotation',0)

subplot(3,1,3)
plot(xvec,tau(3,:),'k','linewidth',1)
hold on
plot(SPess(1,:)',0*SPess(1,:)','k','linewidth',4)
plot(E((RES'<0.0001)&(E>SPess(1,2)+0.001)),E((RES'<0.0001)&(E>SPess(1,2)+0.001))*0,'ok','linewidth',1,'markersize',8)
axis([-2.1,-1.6,0,0.15])

ax=gca; ax.FontSize=16;
xlabel('$x$','interpreter','latex','fontsize',18)
legend({'$\tau_{n_2,1000}(H+V-xI)$','$\mathrm{Sp}_{\mathrm{ess}}(H+V)$','$\mathrm{Sp}_{\mathrm{d}}(H+V)$'},'interpreter','latex','fontsize',14)
text(-2.05,0.1,'$n_2=2$','interpreter','latex','fontsize',18,'Rotation',0)
ff.Position=[360.0000   97.6667  560.0000*1.2  420.0000*1.1];

%% Output of algorithms for discrete and essential spectra (Figure 8.4)

n2 = 100;
xvec = -3:0.001:4;

tau = zeros(1,length(xvec)); dist = tau;
pf = parfor_progress(length(xvec));
pfcleanup = onCleanup(@() delete(pf));
warning('off','all')
for jj=1:length(xvec)
    Indx = find(abs(nn)>n2-1);
    tau(jj) = svds(H(:,Indx)-xvec(jj)*I(:,Indx),1,'smallest');
    dist(jj) = svds(H-xvec(jj)*I,1,'smallest');
    parfor_progress(pf);
end
warning('on','all')


%% Plot the results

ff=figure;

subplot(3,1,1)
sigma_ess = xvec(tau<=1/n2);
stem(sigma_ess,0*sigma_ess+1,'k','markersize',0.0001,'linewidth',0.5)
ax=gca; ax.FontSize=16;
xlim([-3,4])
ylim([0,1])
box on
set(gca,'ytick',[])
title('Algorithm for Essential Spectrum','interpreter','latex','fontsize',16) 

subplot(3,1,2)
sigma_d = xvec((dist<=1/(n1))&(tau>1/n2-2/(n1)));
stem(sigma_d,0*sigma_d+1,'k','markersize',0.0001,'linewidth',0.5)
xlim([-3,4])
ff.Position=[360.0000   97.6667  560.0000*1.2  420.0000];
ax=gca; ax.FontSize=16;
xlim([-3,4])
set(gca,'ytick',[])
title('Algorithm for Discrete Spectrum','interpreter','latex','fontsize',16)

subplot(3,1,3)
stem(diag(D),0*diag(D)+1,'k','markersize',0.0001,'linewidth',0.5)
ax=gca; ax.FontSize=16;
xlim([-3,4])
set(gca,'ytick',[])
title('Finite Section','interpreter','latex','fontsize',16)
text(-2.3,1.2,'spurious','interpreter','latex','fontsize',14,'Rotation',0)
annotation('arrow',[0.29 0.315],[0.19 0.1]+0.1,'linewidth',2)


%% Compute multiplicity (Figure 8.5)

E=diag(D);
xvec2=[-2.2; E(abs(E+1.3)==min(abs(E+1.3)));E(abs(E+1.7)==min(abs(E+1.7)))]; xvec2(2)=[];
n2vec2=0:9;
h = zeros(length(n2vec2),length(xvec2));
warning('off','all')
for ii=1:length(n2vec2)
    for jj=1:length(xvec2)
        h(ii,jj) = sum(max(1-(n2vec2(ii)+1)*svds(H-xvec2(jj)*I,n2vec2(ii)+2,'smallest'),0));
    end
end
warning('on','all')


%% Plot the results

figure

plot(n2vec2+1,round(h(:,1)),'.-k','markersize',20,'linewidth',2);
hold on
plot(n2vec2+1,(h(:,3)),'.--k','markersize',20,'linewidth',1.5);
plot(n2vec2+1,h(:,2),'.:k','markersize',20,'linewidth',2);
legend({'$x\approx -2.2000$','$x\approx -1.6892$','$x\approx -1.3018$'},'interpreter','latex','fontsize',16,'location','best')
title('$h_{n_2,1000}(x,A)$','interpreter','latex','fontsize',18)
xlabel('$n_2$','interpreter','latex')
ax=gca; ax.FontSize=18;
xlim([1,10])

ff=figure;
subplot(3,1,2)
stem(sigma_ess,0*sigma_ess+1,'k','markersize',0.0001,'linewidth',0.5)
hold on
stem(sigma_d,0*sigma_d+1,'k','markersize',0.0001,'linewidth',0.5)
ax=gca; ax.FontSize=16;
xlim([-3,4])
ylim([0,1])
box on
set(gca,'ytick',[])
title('Spectrum','interpreter','latex','fontsize',16) 

subplot(3,1,1)
xlim([-3,4])
ylim([0,1])
annotation('arrow',1.53/7*[1 1],[0.1 0]+0.582,'linewidth',2)
hold on
annotation('arrow',1.92/7*[1 1],[0.1 0]+0.582,'linewidth',2)
annotation('arrow',2.22/7*[1 1],[0.1 0]+0.582,'linewidth',2)
box on
set(gca,'ytick',[])
ax=gca; ax.FontSize=16;
xlim([-3,4])
ylim([0,1])
title('Spectrum','interpreter','latex','fontsize',16) 

subplot(3,1,3)
xlim([-3,4])
ylim([0,1])
box on
set(gca,'ytick',[])
ax=gca; ax.FontSize=16;
xlim([-3,4])
ylim([0,1])
title('Spectrum','interpreter','latex','fontsize',16) 
ff.Position=[360.0000   97.6667  560.0000*1.2  420.0000];



