clear
close all

% Add utils and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'utils');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'data_online');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%% CODE FOR INTEGRAL OPERATOR EXAMPLE %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Adaptive selection of discretisation sizes
X = -2.5:0.001:2.5;                      %evaluation pts
f = @(x) (exp(x)+x)/sqrt(5.765044839199454);                         %measure wrt f(x)
a = {@(x) x, @(x,y) cos(2*x.*y)};       %integral op. coeffs

mu1=intMeas(a,f,X,0.1,'Order',1);           %epsilon = 0.1
mu2=intMeas(a,f,X,0.01,'Order',1);          %epsilon = 0.01
mu3=intMeas(a,f,X,0.001,'Order',1);         %epsilon = 0.001

figure
semilogy(X,mu1,'k:','linewidth',2)
hold on
semilogy(X,mu2,'k--','linewidth',2)
semilogy(X,mu3,'k','linewidth',2)
xlim([X(1) X(end)])
ylim([10^(-5) 10^3])
legend({'$\epsilon=0.1$','$\epsilon=0.01$','$\epsilon=0.001$'},'fontsize',16,'interpreter','latex','location','north')

xlabel('$x$','interpreter','latex','fontsize',18)
title('$\mu_v^\epsilon(x)$','interpreter','latex','fontsize',18)
yticks(10.^(-4:2:2))
ax=gca; ax.FontSize=18;

%% Fixed discretisation sizes
mu1=intMeas(a,f,X,0.1,'Order',1,'N',200);           %epsilon = 0.1
mu2=intMeas(a,f,X,0.01,'Order',1,'N',200);          %epsilon = 0.01
mu3=intMeas(a,f,X,0.001,'Order',1,'N',200);         %epsilon = 0.001

figure
semilogy(X,mu1,'k:','linewidth',2)
hold on
semilogy(X,mu2,'k--','linewidth',2)
semilogy(X,mu3,'k','linewidth',2)
xlim([X(1) X(end)])
ylim([10^(-5) 10^3])
ax = gca; ax.FontSize = 14;
xlabel('$x$','interpreter','latex','fontsize',18)
title('Approximation with fixed $N=200$','interpreter','latex','fontsize',18)
yticks(10.^(-4:2:2))
ax=gca; ax.FontSize=18;


a2 = axes();
a2.Position = [0.42 0.73 0.3 0.15]; % xlocation, ylocation, xsize, ysize
semilogy(X,mu3,'k','LineWidth',1.5)
hold on
semilogy(X,mu2,'k--','LineWidth',1.5)
semilogy(X,mu1,'k:','LineWidth',1.5)
axis([-0.5,-0.3,0.005,0.5])

%% Interior layers
[u1,ccn1] = intMeas3(a,f,0.5,0.1,'Order',1);
[u2,ccn2] = intMeas3(a,f,0.5,0.05,'Order',1);
[u3,ccn3] = intMeas3(a,f,0.5,0.01,'Order',1);

figure
plot(ccn3,imag(u3),'g','LineWidth',1)
hold on
plot(ccn2,imag(u2),'r','LineWidth',1.5)
plot(ccn1,imag(u1),'k','LineWidth',2)
legend({'$\epsilon=0.01$','$\epsilon=0.05$','$\epsilon=0.1$'},'fontsize',16,'interpreter','latex','location','southwest')
xlabel('$x$','interpreter','latex','fontsize',18)
title('$u^\epsilon(x)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;


%% Errors in terms of N
Nvec=100:100:4000;
U_ref1=intMeas(a,f,0.5,0.05,'Order',1,'N',10^4);
U_ref2=intMeas(a,f,0.5,0.01,'Order',1,'N',10^4);
U_ref3=intMeas(a,f,0.5,0.005,'Order',1,'N',10^4);
E=zeros(length(Nvec),3);

for j=1:length(Nvec)
    
    u1=intMeas(a,f,0.5,0.05,'Order',1,'N',Nvec(j));
    u2=intMeas(a,f,0.5,0.01,'Order',1,'N',Nvec(j));
    u3=intMeas(a,f,0.5,0.005,'Order',1,'N',Nvec(j));
    
    E(j,1)=abs(u1-U_ref1)/abs(U_ref1);
    E(j,2)=abs(u2-U_ref2)/abs(U_ref2);
    E(j,3)=abs(u3-U_ref3)/abs(U_ref3);
end


figure
semilogy(Nvec,E(:,1),'marker','o','MarkerSize',8,'color','k')
hold on
semilogy(Nvec,E(:,2),'marker','x','MarkerSize',8,'color','k')
semilogy(Nvec,E(:,3),'marker','.','MarkerSize',20,'color','k')
axis([0,3000,10^(-15),1])

text(1500,10^(-6),'$\epsilon=0.005$','interpreter','latex','fontsize',16,'Rotation',-35)
text(800,10^(-7),'$\epsilon=0.01$','interpreter','latex','fontsize',16,'Rotation',-50)
text(300,10^(-8),'$\epsilon=0.05$','interpreter','latex','fontsize',16,'Rotation',-78)

xlabel('$N$','interpreter','latex','fontsize',18)
title('$| \mu_{v,N}^\epsilon(x_0)-\mu_v^\epsilon(x_0)|/|\mu_v^\epsilon(x_0)|$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;


%% Errors in terms of epsilon

ep_vec=10.^(-2:0.1:0);
U_ref1 = intMeas(a,f,0.5,0.005,'Order',6);
E = zeros(length(ep_vec),6);

for j=1:length(ep_vec)
   u=zeros(6,1);
   for jj=1:6
        u=intMeas(a,f,0.5,ep_vec(j),'Order',jj);
        E(j,jj)=abs(u-U_ref1)/abs(U_ref1);
   end
end


figure
loglog(ep_vec,E(:,1),'marker','o','MarkerSize',8,'color','k','linewidth',1.5)
hold on
e2=ep_vec(3:10)/1000;
loglog(ep_vec(3:10),e2.*log(1+1./e2)*200,'--k','linewidth',2)
ylim([0.01,1])
text(0.016,0.06,'$\mathcal{O}(\epsilon\log(1+\epsilon^{-1}))$','interpreter','latex','fontsize',16,'Rotation',33)
xlabel('$\epsilon$','interpreter','latex','fontsize',18)
title('$|\rho_v(x_0)-\mu_v^\epsilon(x_0)|/|\rho_v(x_0)|$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

figure
loglog(ep_vec,E(:,1),'marker','o','MarkerSize',8,'color','k','linewidth',1.5)
hold on
loglog(ep_vec,E(:,2),'marker','s','MarkerSize',8,'color','k','linewidth',1.5)
loglog(ep_vec,E(:,3),'marker','^','MarkerSize',8,'color','k','linewidth',1.5)
loglog(ep_vec,E(:,4),'marker','.','MarkerSize',25,'color','k','linewidth',1.5)
loglog(ep_vec,E(:,5),'marker','*','MarkerSize',8,'color','k','linewidth',1.5)
loglog(ep_vec,E(:,6),'marker','+','MarkerSize',8,'color','k','linewidth',1.5)
legend({'$m=1$','$m=2$','$m=3$','$m=4$','$m=5$','$m=6$'},'fontsize',16,'interpreter','latex','location','southeast')
xlabel('$\epsilon$','interpreter','latex','fontsize',18)
title('$|\rho_v(x_0)-[K_\epsilon\ast\mu_v](x_0)|/|\rho_v(x_0)|$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

%% Zoomed in plot
X=0.8:0.001:1.2;
mu1=intMeas(a,f,X,0.01,'Order',1);           %epsilon = 0.1
mu2=intMeas(a,f,X,0.01,'Order',2);          %epsilon = 0.01
mu3=intMeas(a,f,X,0.01,'Order',6);         %epsilon = 0.001

figure
semilogy(X,mu1,'k:','linewidth',2)
hold on
semilogy(X,mu2,'k--','linewidth',2)
semilogy(X,abs(mu3),'k','linewidth',2)
legend({'$m=1$','$m=2$','$m=6$'},'fontsize',16,'interpreter','latex','location','southwest')
xlim([X(1) X(end)])
ylim([10^(-12) 10^0])
xlabel('$x$','interpreter','latex','fontsize',18)
title('$\left|[K_\epsilon\ast\mu_v](x)\right|$','interpreter','latex','fontsize',18)
yticks(10.^(-15:5:2))
ax=gca; ax.FontSize=18;

