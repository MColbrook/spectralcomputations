close all
clear

%% Example densities

f = @(r) cos(r).*exp(-r.^2/10)/sqrt(0.997507996719327);
X=-1:0.05:4;
epsilon = 0.2;

V={@(r) 0*r, @(r) exp(-r)-1/2, 1};
mu1=rseMeas(V,f,X,epsilon,'Order',6);
V={@(r) 0*r, @(r) exp(-r)-1/2, 2};
mu2=rseMeas(V,f,X,epsilon,'Order',6); 
V={@(r) 0*r, @(r) exp(-r)-1/2, 3};
mu3=rseMeas(V,f,X,epsilon,'Order',6); 


figure
plot(X,max(mu1,10^(-16)),'k','linewidth',2)
hold on
plot(X,max(mu2,10^(-16)),'k--','linewidth',2)
plot(X,max(mu3,10^(-16)),'k:','linewidth',2)
ylim([0,0.17])

legend({'$\ell=1$','$\ell=2$','$\ell=3$'},'fontsize',16,'interpreter','latex','location','northeast')

xlabel('$x$','interpreter','latex','fontsize',18)
title('$[K_\epsilon\ast\mu_v](x)$','interpreter','latex','fontsize',18)
ax = gca; ax.FontSize = 18;

%% Compute the probability using quadrature
V={@(r) 0*r, @(r) exp(-r)-1/2, 1};
Nvec = [1:20,30];

for N = Nvec
    [ccn,ccw]=chebpts(N,[1 3]);
    mu=rseMeas(V,f,ccn,0.1,'Order',6);
    sum(mu(:).*ccw(:))
end
