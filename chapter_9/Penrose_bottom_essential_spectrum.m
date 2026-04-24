clear
close all

% Add utils and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'utils');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'data_online');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%%%%%%%%% CODE FOR FIGURE 9.6 %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Build the matrix of the operator

k = 14; % level of approximation
TYPE = 3; % tile type

% build the matrix

if TYPE == 1
    [L, ~, loc, ~, indx_bdy] = circ_lap_tiling(k);
    X = 0:0.01:7;
elseif TYPE == 2
    [L, loc, ~, indx_bdy] = circ_lap_dual(k);
    X = 0:0.01:13;
else
    [L, loc, ~, indx_bdy] = circ_lap_dual_rhomb(k);
    X = 0:0.01:10;
end

[~,I]=sort_mi(loc);
Id = speye(length(I));
pt(I) = 1:length(I);
indx_bdy = pt(indx_bdy);
L = L(I,I);
L(:,indx_bdy)=[];
Id(:,indx_bdy)=[]; 


%% Approximate the bottom of the essential spectrum

nvec = round(10.^(2:1:6));
ct = 1;
ess = zeros(100,5);

for N = nvec % see output as it is produced

    L2 = L(:,1:N);
    Id2 = Id(:,1:N);
    I = find(sum(abs(L2)+abs(Id2),2)<0.1);
    L2(I,:) = [];
    Id2(I,:) = [];
    
    ess(:,ct) = sort(svds(L2+Id2,100,'smallest'),'ascend')-1;
    semilogy(ess)
    xlabel('$n_2$','interpreter','latex','fontsize',18)
    title('$\sigma_{\mathrm{inf}+n_2-1}(\mathcal{P}_{f(n_1)}L_{T_3}\mathcal{P}^*_{n_1})$','interpreter','latex','fontsize',18)
    ct = ct+1;
    pause(1)
end

%% Plot the results

close all
figure
semilogy(ess(:,1),'k','linewidth',3)
hold on
plot(ess(:,2),':k','linewidth',2)
plot(ess(:,3),'-.k','linewidth',2)
plot(ess(:,4),'--k','linewidth',2)
plot(ess(:,5),'-k','linewidth',1)
xlabel('$n_2$','interpreter','latex','fontsize',18)
    title('$\sigma_{\mathrm{inf}+n_2-1}(\mathcal{P}_{f(n_1)}L_{T_3}\mathcal{P}^*_{n_1})$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=17;
ylim([10^(-5),10])
legend({'$n_1=10^2$','$n_1=10^3$','$n_1=10^4$','$n_1=10^5$','$$n_1=10^6$'},'interpreter','latex',...
    'fontsize',16,'location','southeast','NumColumns',2)









