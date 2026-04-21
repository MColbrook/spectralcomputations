clear
close all

%% Plot basis and condition numbers (uses Chebfun)

nvec = 0.25:0.25:10;
cc = 0*nvec;
N = 100;

ct =1;
for n = nvec
    V = chebfun(@(x) exp(-x).*sin((1:N)*pi/2*(x/n+1)),[-n,n],'turbo');
    for jj = 1:N
        V(:,jj) = V(:,jj)/norm(V(:,jj));
    end
    cc(ct)=cond(V);
    ct = ct+1;
end

figure
semilogy(nvec,cc,'.k','markersize',15)
xlabel('$n$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
exportgraphics(gcf,'CHAP1_bad_basis2.pdf','ContentType','vector','BackgroundColor','none','Colorspace','gray')

figure
V = chebfun(@(x) exp(-x).*sin((1:10)*pi/2*(x/5+1)),[-5,5],'turbo');
plot(V,'k','linewidth',1)
xlabel('$x$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
grid on
exportgraphics(gcf,'CHAP1_bad_basis1.pdf','ContentType','vector','BackgroundColor','none','Colorspace','gray')

%% Plot different spectra

k=-100:0.01:100;

n=2;

E2 = 1:1000;
E2 = 1+E2.^2*pi^2/(4*n^2);

plot(real(E2),imag(E2),'k.','markersize',15)
hold on
plot(k.^2,2*k,'k','linewidth',1)
xlim([0,100])
ylim([-22,22])
text(40,-5,'$n=2$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
ax = gca; ax.FontSize = 18;

exportgraphics(gcf,'CHAP1_davies1.pdf','ContentType','vector','BackgroundColor','none','Colorspace','gray')

n=10;

E2 = 1:1000;
E2 = 1+E2.^2*pi^2/(4*n^2);
figure
plot(real(E2),imag(E2),'k.','markersize',15)
hold on
plot(k.^2,2*k,'k','linewidth',1)
xlim([0,100])
ylim([-22,22])
text(40,-5,'$n=10$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
ax = gca; ax.FontSize = 18;

exportgraphics(gcf,'CHAP1_davies2.pdf','ContentType','vector','BackgroundColor','none','Colorspace','gray')



