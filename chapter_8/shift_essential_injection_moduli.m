clear
close all

%%%%%%%%%%%%%%%%%%% FIGURE 8.2 (SHIFT EXAMPLE) %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Standard injection moduli

xvec= 0:0.005:2; % grid of points over which we approximate injection moduli

n1 = 100;
n2 = 0;
sigma = 0*xvec;

A = [sparse(zeros(n1,1)),speye(n1,n1);sparse(zeros(1,n1+1))];
I = speye(n1,n1);

A = A(1:n1,(n2+1):n1); I = I(1:n1,(n2+1):n1);

for j = 1:length(xvec)
    sigma(j) = svds(A-xvec(j)*I,1,'smallest');
end

figure
tiledlayout(3,1,"TileSpacing","compact");
nexttile([2,1])
plot(xvec,sigma,'k','linewidth',2)
xlim([xvec(1),xvec(end)])
ylim([0,1])
set(gca,'xtick',[])
title('$\sigma_{\mathrm{inf}}((A-zI)\mathcal{P}_{100}^*)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
nexttile
stem(0:0.01:1,0*(0:0.01:1)+1,'k','markersize',0.0001,'linewidth',2)
title('$\mathrm{Sp}(A)$','interpreter','latex','fontsize',18)
xlim([xvec(1),xvec(end)])
ylim([0.25,0.75])
ax=gca; ax.FontSize=18;
set(gca,'ytick',[])
xlabel('$|z|$','interpreter','latex','fontsize',18)


%% Essential injection moduli (with the Q projection)

clear
xvec= 0:0.005:2;

n1 = 100;
n2 = 1;
tau = 0*xvec;

A = [sparse(zeros(n1,1)),speye(n1,n1);sparse(zeros(1,n1+1))]; %A = A';
I = speye(n1,n1);

A = A(1:n1,(n2+1):n1); I = I(1:n1,(n2+1):n1);

for j = 1:length(xvec)
    tau(j) = svds(A-xvec(j)*I,1,'smallest');
end


figure
tiledlayout(3,1,"TileSpacing","compact");
nexttile([2,1])
plot(xvec,tau,'k','linewidth',2)
xlim([xvec(1),xvec(end)])
ylim([0,1])
set(gca,'xtick',[])
title('$\sigma_{\mathrm{inf}}((A-zI)\mathcal{Q}_{1,100}^*)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
nexttile
stem(1,1,'k','markersize',0.0001,'linewidth',2)
title('$\mathrm{Sp}_{\mathrm{ess},2}(A)$','interpreter','latex','fontsize',18)
xlim([xvec(1),xvec(end)])
ax=gca; ax.FontSize=18;
set(gca,'ytick',[])
ylim([0.25,0.75])
xlabel('$|z|$','interpreter','latex','fontsize',18)


