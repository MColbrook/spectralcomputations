clear
close all

% Add Koopman algorithms and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'Koopman_algorithms');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'Koopman_datasets');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%%%%% CODE FOR TURBULENCE EXAMPLE %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data set (see link in README file)
load('turbulent_data')
LL = 3000; % number of timesteps for prediction, increase to 5000 for book plot
sam_num = 2; % number of random initial conditions for TKE plot, increase to 100 for book plot

%% Compute measurement matrices
[~,~,V] = svd(transpose(DATA(:,ind1))/sqrt(M),'econ');

PX = transpose(DATA(:,ind1))*V(:,1:N);
PY = transpose(DATA(:,ind1+1))*V(:,1:N);

% mpEDMD
[mpK,mpV,mpD] = mpEDMDqr(PX,PY,1/M);
mpLAM = diag(mpD);

% EDMD
K = PX\PY;
[V,LAM]=eig(K,'vector');

%% Koopman mode decomposition
sigu = DATA(1:51150,ind1);
sigv = DATA(51151:end,ind1);

kmd = single([(PX*V)\transpose(sigu),(PX*V)\transpose(sigv)]);
mpkmd = single([(PX*mpV)\transpose(sigu),(PX*mpV)\transpose(sigv)]);

ri = randsample(length(ind1),sam_num); % random initial conditions to sample over
LAM = single(LAM);
mpLAM = single(mpLAM);
PX = single(PX);
V = single(V);
mpV = single(mpV);

%% Turbulent kinetic energy
kmd1 = kmd(:,1:51150);
kmd2 = kmd(:,51151:end);

mpkmd1 = mpkmd(:,1:51150);
mpkmd2 = mpkmd(:,51151:end);

pf = parfor_progress(LL*length(ri));
pfcleanup = onCleanup(@() delete(pf));

for jjj=1:length(ri)
    X = PX(ri(jjj),:)*V;
    mpX = PX(ri(jjj),:)*mpV;
    
    uf = single(zeros(310,165,LL));
    vf = uf;
    mpuf = uf;
    mpvf = uf;

    for j=1:LL
        t = real((X.*transpose(LAM(:).^(j-1)))*kmd1);
        uf(:,:,j) = reshape(t,[310,165]);
        t = real((X.*transpose(LAM(:).^(j-1)))*kmd2);
        vf(:,:,j) = reshape(t,[310,165]);

        t = real((mpX.*transpose(mpLAM(:).^(j-1)))*mpkmd1);
        mpuf(:,:,j) = reshape(t,[310,165]);
        t = real((mpX.*transpose(mpLAM(:).^(j-1)))*mpkmd2);
        mpvf(:,:,j) = reshape(t,[310,165]);
        parfor_progress(pf);
    end

    TKE_EDMD0 = flow_TKE(uf,vf);
    TKE_mpEDMD0 = flow_TKE(mpuf,mpvf);

    if jjj==1
        TKE_EDMD=TKE_EDMD0; TKE_mpEDMD=TKE_mpEDMD0;           
    else
        TKE_EDMD=(TKE_EDMD*(jjj-1)+TKE_EDMD0)/jjj; TKE_mpEDMD=(TKE_mpEDMD*(jjj-1)+TKE_mpEDMD0)/jjj; 
    end

end


%% TKE plots

III=find(abs(y-5)==min(abs(y-5)));
tt = squeeze(TKE_mpEDMD(:,III,:)); tt = mean(tt(:));

figure
semilogy((1:LL)/1000,squeeze(mean(mean(TKE_mpEDMD(:,III,:),2),1))/2500,'color',[1,1,1]*0.6,'linewidth',1);
hold on
plot((1:LL)/1000,squeeze(mean(mean(TKE_EDMD(:,III,:),2),1))/2500,'-k','linewidth',2)
semilogy([0,LL/1000],tt/2500 +zeros(2,1),'--k','linewidth',2)
legend({'mpEDMD','EDMD','Flow'},'interpreter','latex','fontsize',20,'location','northwest')
xlabel('Time (s)','interpreter','latex','fontsize',20)
title('$y\approx 5$mm','interpreter','latex','fontsize',20)
ylim([0.001/2,1])
xlim([0,5])
ax = gca; ax.FontSize = 20;


III=find(abs(y-35)==min(abs(y-35)));
tt = squeeze(TKE_mpEDMD(:,III,:)); tt = mean(tt(:)); 

figure
semilogy((1:LL)/1000,squeeze(mean(mean(TKE_mpEDMD(:,III,:),2),1))/2500,'color',[1,1,1]*0.6,'linewidth',1);
hold on
plot((1:LL)/1000,squeeze(mean(mean(TKE_EDMD(:,III,:),2),1))/2500,'-k','linewidth',2)
semilogy([0,LL/1000],tt/2500 +zeros(2,1),'--k','linewidth',2)
legend({'mpEDMD','EDMD','Flow'},'interpreter','latex','fontsize',20,'location','northwest')
xlabel('Time (s)','interpreter','latex','fontsize',20)
title('$y\approx 35$mm','interpreter','latex','fontsize',20)
ylim([0.001/2,1])
xlim([0,5])
ax = gca; ax.FontSize = 20;


function TKE = flow_TKE(uf,vf)
TKE = (uf.*uf+vf.*vf)/2; % turbulent kinectic energy
end







