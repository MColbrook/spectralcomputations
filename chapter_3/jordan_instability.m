%% Jordan matrix pseudospectra

close all
clear

x_pts=-1:0.02:1;    y_pts=-1:0.02:1;
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);		% complex points where we compute pseudospectra
RES=0*z_pts+10;

pf = parfor_progress(length(z_pts));
pfcleanup = onCleanup(@() delete(pf));

k = 5;
J = eye(k+1); J(:,end)=[]; J(1,:)=[];
for jj=1:length(z_pts)
    B = J - z_pts(jj)*eye(size(J));
    RES(jj) = min(svd(B));
    parfor_progress(pf);
end
RES=reshape(RES,length(y_pts),length(x_pts));

%% Plot pseudospectra

v=(10.^(-4:0.4:0));

figure
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(max(real(RES),min(v))),log10(v),'LineColor','k',...
    'linewidth',0.5,'ShowText','off');%,'FaceAlpha',0.6);
cbh=colorbar;
cbh.Ticks=log10(10.^(-4:1:0));
cbh.TickLabels=10.^(-4:1:0);
clim([-4,0]);
colormap gray
axis equal;

title('$\mathrm{Sp}_{\epsilon}(J_5)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
box on

%% Instability plot

x2 = [0.05,0.2,0.5,0.9];
RES2 = zeros(length(x2),10);

for ii = 1:10
    J = eye(ii+1); J(:,end)=[]; J(1,:)=[];
    for jj=1:length(x2)
        B = J - x2(jj)*eye(size(J));
        RES2(jj,ii) = min(svd(B));
    end
end

figure
semilogy(1:10,1./RES2,'k','linewidth',2)
title('$\|(J_k-zI)^{-1}\|$','interpreter','latex','fontsize',18)
xlabel('$k$','interpreter','latex','fontsize',18)
box on
ax=gca; ax.FontSize=18;
xlim([1,10])
text(7,10./RES2(1,7),'$|z|=0.05$','interpreter','latex','fontsize',16,'Rotation',30)
text(7,10./RES2(2,7),'$|z|=0.2$','interpreter','latex','fontsize',16,'Rotation',18)
text(7,5./RES2(3,7),'$|z|=0.5$','interpreter','latex','fontsize',16,'Rotation',10)
text(7,5./RES2(4,5),'$|z|=0.9$','interpreter','latex','fontsize',16,'Rotation',3)
exportgraphics(gcf, 'CHAP3_Jordan2.pdf','ContentType','vector','BackgroundColor','none','Resolution',300,'Colorspace','gray')

