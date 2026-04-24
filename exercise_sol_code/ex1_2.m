% Hermite function discretisation of convection–diffusion operator
clear
close all

%% Build large matrix
N = 1000;
bb = spdiags(sqrt((0:N+100)'/2),1,N+100,N+100);
dx = bb-bb';

dx = dx(1:(N+99),1:(N+99)); % derivative operator

L = -dx*dx -2*dx;
L = L(1:N,1:N);

%% Eigenvalues of finite section
Nvec = [10,50,100,1000];
k = -100:0.01:100; spec = k.^2 + 2i*k;

for N = Nvec
    figure
    E = eig(full(L(1:N,1:N)));
    plot(real(E),imag(E),'k.','markersize',15)
    hold on
    plot(real(spec),imag(spec),'k','linewidth',1)
    xlim([0,100])
    ylim([-22,22])

    title(sprintf('$N=%d$',N),'interpreter','latex','fontsize',18)
    xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
    ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
    ax = gca; ax.FontSize = 18;
    exportgraphics(gcf,sprintf('hermite_con_diff%d.pdf',N),'ContentType','vector','BackgroundColor','none','Colorspace','gray')
end



