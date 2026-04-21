clear
close all

%% Set up the matrix
N = 200;
bb=spdiags(sqrt((0:N+100)'/2),1,N+100,N+100);
mx=bb+bb';
dx=bb-bb';

Lap = dx*dx;
H = -Lap + 1i*mx*mx*mx;

%% Eigenvalues and approximate condition numbers: warning condition number plot in book requires larger truncation size and use of higher precision
N = 120;
[V,E] = eig(full(H(1:N,1:N)));
E = diag(E);
[~,I] = sort(real(E)+abs(imag(E))*10^4,'ascend');
E = E(I); V = V(:,I);
cn = transpose(vecnorm(V).^2)./abs(diag(conj(V)'*V));
En = 1:12; cn = cn(1:12);

figure
loglog(En,cn.*exp(-pi/sqrt(3)*(1:En(end))'),'k.','markersize',12)
hold on
loglog(En(5:end),.2./(En(5:end).^(1/4)),'--k','linewidth',2)
text(7,0.135,'$\mathcal{O}(j^{-1/4})$','interpreter','latex','fontsize',18,'Rotation',-26)
legend({'$\|\mathcal{Q}_j\|\exp(-j\pi/\sqrt{3})$'},'fontsize',16,'interpreter','latex','location','northeast')
xlabel('$j$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
exportgraphics(gcf,'CHAP3_HB_cond.pdf','ContentType','vector','BackgroundColor','none','Resolution',300)

%% Compute pseudospectra of finite sections
N = 100;
xpts = 0:1:80;    ypts = 0:1:100;
zpts=kron(xpts,ones(length(ypts),1))+1i*kron(ones(1,length(xpts)),ypts(:));    zpts=zpts(:);		% complex points where we compute pseudospectra
RES=0*zpts;

pf = parfor_progress(length(zpts));
pfcleanup = onCleanup(@() delete(pf));

for jj=1:length(zpts)
    B = H(1:N,1:N) - speye(N,N)*zpts(jj);
    RES(jj) = svds(B,1,'smallest');
    parfor_progress(pf);
end
RES=reshape(RES,length(ypts),length(xpts));

%% Compute pseudospectra using rectangular sections
RES2=0*zpts;

pf = parfor_progress(length(zpts));
pfcleanup = onCleanup(@() delete(pf));
for jj=1:length(zpts)
    B = H(1:(N+4),1:N) - speye(N+4,N)*zpts(jj);
    RES2(jj) = svds(B,1,'smallest');
    parfor_progress(pf);
end
RES2=reshape(RES2,length(ypts),length(xpts));

%% Plot pseudospectra
E = eig(full(H(1:N,1:N)));
v=(10.^(-20:0.5:0));
figure
hold on
contourf(reshape(real(zpts),length(ypts),length(xpts)),reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES),min(v))),log10(v),'LineColor','k',...
    'linewidth',1,'ShowText','off');
hold on
contourf(reshape(real(zpts),length(ypts),length(xpts)),-reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES),min(v))),log10(v),'LineColor','k',...
    'linewidth',1,'ShowText','off');
cbh=colorbar;
cbh.Ticks=log10(10.^(-16:2:0));
cbh.TickLabels=["1e-16","1e-14","1e-12","1e-10","1e-08","1e-06","1e-04","1e-02","1"];
clim([-16,0])
colormap gray
ax=gca; ax.FontSize=14;
axis tight
plot(real(E),imag(E),'k.','markersize',5)
axis([min(xpts),max(xpts),min(ypts),max(ypts)])
xlim([0,80])
ylim([-100,100])
title('Finite Section','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
box on


v=(10.^(-20:0.5:0));
figure
hold on
contourf(reshape(real(zpts),length(ypts),length(xpts)),reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES2),min(v))),log10(v),'LineColor','k',...
    'linewidth',1,'ShowText','off');
contourf(reshape(real(zpts),length(ypts),length(xpts)),-reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES2),min(v))),log10(v),'LineColor','k',...
    'linewidth',1,'ShowText','off');
cbh=colorbar;
cbh.Ticks=log10(10.^(-16:2:0));
cbh.TickLabels=["1e-16","1e-14","1e-12","1e-10","1e-08","1e-06","1e-04","1e-02","1"];
clim([-16,0])
colormap gray
ax=gca; ax.FontSize=14;
axis tight
hold on
axis([min(xpts),max(xpts),min(ypts),max(ypts)])
xlim([0,80])
ylim([-100,100])
title('\texttt{PseudoSpec}','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)

ax=gca; ax.FontSize=18;
box on


%% Distance plots
xpts=0.8:0.001:1.5;
L2=length(xpts);
dist1=xpts*0;
dist2=dist1;
dist3=dist1;

g = (exp(pi/sqrt(3)*(1:100))-1)/(exp(pi/sqrt(3))-1).*exp(3*(1:100).^(6/5));

for j=1:L2
    B=H(:,1:N)-xpts(j)*speye(size(H,1),N);
    dist1(j)=svds(B,1,'smallest');
    dist2(j)=min(abs(E-xpts(j)));
    dist3(j)=dist1(j)*1*cn(1);
end

figure
plot(xpts,dist2,'k','linewidth',3)
hold on
loglog(xpts,dist1,'--k','linewidth',2)
plot(xpts,dist3,':k','linewidth',2)
xlim([xpts(1),xpts(end)])
legend({'$\mathrm{dist}(x,\mathrm{Sp}(H))$','$\gamma(x,H)$','$\|\mathcal{Q}_1\|\gamma(x,H)$'},...
    'fontsize',16,'interpreter','latex','location','north')
xlabel('$x$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

exportgraphics(gcf,'CHAP3_HB_gamma1.pdf','ContentType','vector','BackgroundColor','none','Resolution',300,'Colorspace','gray')

figure
plot(xpts,dist1-dist2,'--k','linewidth',2)
hold on
plot(xpts,dist3-dist2,':k','linewidth',2)
ax=gca; ax.FontSize=14;
xlim([xpts(1),xpts(end)])
legend({'$\gamma(x,H)-\mathrm{dist}(x,\mathrm{Sp}(H))$','$\|\mathcal{Q}_1\|\gamma(x,H)-\mathrm{dist}(x,\mathrm{Sp}(H))$'},'fontsize',16,'interpreter','latex','location','south')
xlabel('$x$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

exportgraphics(gcf,'CHAP3_HB_gamma2.pdf','ContentType','vector','BackgroundColor','none','Resolution',300,'Colorspace','gray')

%% Plot eigenfunction, warning larger j requires larger truncation size and use of higher precision
XX = -5:0.005:5; XX = XX(:);
HH = zeros(length(XX),size(V,1));
HH(:,1)=exp(-XX.^2/2)/(pi^(1/4));
HH(:,2)=exp(-XX.^2/2)/(pi^(1/4)).*XX*sqrt(2);

for j = 3:size(V,1)
    n = j-1;
    HH(:,j)=HH(:,j-1).*XX*sqrt(2/n)-HH(:,j-2)*sqrt((n-1)/n);
end


figure
Vx = HH*V(:,10);
plot(XX,real(Vx),'k','linewidth',2)
ax=gca; ax.FontSize=20;
xlabel('$x$','interpreter','latex','fontsize',24)
title('$\mathrm{Re}(u_j(x)),j=10$','interpreter','latex','fontsize',24)

figure
plot(XX,real(Vx.^2),'k','linewidth',2)
ax=gca; ax.FontSize=20;
xlabel('$x$','interpreter','latex','fontsize',24)
title('$\mathrm{Re}(u_j^2(x)),j=10$','interpreter','latex','fontsize',24)

