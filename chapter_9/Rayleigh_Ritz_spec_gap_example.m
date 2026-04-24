clear
close all

%%%%%%%%%%% CODE FOR SPECTRAL GAP EXAMPLE (FIGURE 9.10) %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Set up the matrix
N=200;
bb=spdiags(sqrt((0:N+100)'/2),1,N+100,N+100);
mx=bb+bb';
dx=bb-bb';

dx=dx(1:(N+99),1:(N+99)); % derivative operator
mx=mx(1:(N+99),1:(N+99)); % multiplication by x

Lap=dx*dx; % Laplacian
T=-Lap+mx*mx*mx*mx; % add quartic potential


%% Convergence of spectral gaps
E=zeros(2,5,100)+NaN;

for jj = 1:100
    E(1,1:min(jj,5),jj) = eigs(T(1:jj,1:jj),min(jj,5),0);
    E(2,1:min(jj,5),jj) = eigs(-Lap(1:jj,1:jj),min(jj,5),0);
end


%% Plot results

figure
loglog(1:100,squeeze(E(1,2,:)-E(1,1,:)),'k','linewidth',3)
hold on
plot(1:100,squeeze(E(1,3,:)-E(1,1,:)),':k','linewidth',2)
plot(1:100,squeeze(E(1,4,:)-E(1,1,:)),'-.k','linewidth',2)
plot(1:100,squeeze(E(1,5,:)-E(1,1,:)),'k','linewidth',1)
ylim([2,45])
xlabel('$n$','interpreter','latex','fontsize',18)
title('$V(x)=x^4$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

figure
loglog(1:100,squeeze(E(2,2,:)-E(2,1,:)),'k','linewidth',3)
hold on
plot(1:100,squeeze(E(2,3,:)-E(2,1,:)),':k','linewidth',2)
plot(1:100,squeeze(E(2,4,:)-E(2,1,:)),'-.k','linewidth',2)
plot(1:100,squeeze(E(2,5,:)-E(2,1,:)),'k','linewidth',1)
xlabel('$n$','interpreter','latex','fontsize',18)
title('$V(x)=0$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
legend({'$\lambda_{2}^{(n)}-\lambda_{1}^{(n)}$','$\lambda_{3}^{(n)}-\lambda_{1}^{(n)}$',...
    '$\lambda_{4}^{(n)}-\lambda_{1}^{(n)}$','$\lambda_{5}^{(n)}-\lambda_{1}^{(n)}$'},'interpreter','latex',...
    'fontsize',16,'location','southwest')


