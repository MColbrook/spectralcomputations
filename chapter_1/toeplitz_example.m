clear
close all

%% Symbol and spectra

a = @(t) -t.^(-4) -(3+2i)*t.^(-3) +1i*t.^(-2) + 1./t +10*t+(3+1i)*t.^2 + 4*t.^3 +1i*t.^4;
z = a(exp(1i*(0:0.001:2*pi)));

%% Finite section and circulant eigenvalues

nn = 50; % truncation size - adjust accordingly

A = sptoeplitz([0;10;3+1i;4;1i;sparse(nn-5,1)],[0,1,1i,-(3+2i),-1,sparse(1,nn-5)]); % finite section
B = sptoeplitz([0;10;3+1i;4;1i;sparse(nn-9,1);-1;-(3+2i);1i;1],[0,1,1i,-(3+2i),-1,sparse(1,nn-9),1i,4,3+1i,10]); % circulant

figure % finite section
plot(real(z),imag(z),'k','linewidth',1)
E = eig((full((A)))); E=sort(E);
hold on
plot(real(E),imag(E),'k.','markersize',15)
axis([-11,17,-16,17])
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
ax = gca; ax.FontSize = 18;

figure % circulant
plot(real(z),imag(z),'k','linewidth',1)
E2 = eig((full((B)))); E2=sort(E2);
hold on
plot(real(E2),imag(E2),'k.','markersize',15)
axis([-11,17,-16,17])
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
ax = gca; ax.FontSize = 18;

%% Compute pseudospectra

x_pts=-11:0.2:17;    y_pts=-16:0.2:17;
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);		% complex points where we compute pseudospectra
RES = 0*z_pts + 10;
RES2 = RES;

pf = parfor_progress(length(z_pts));
pfcleanup = onCleanup(@() delete(pf));

for jj=1:length(z_pts)
    RES(jj) = svds(A - z_pts(jj)*speye(size(A)),1,'smallest');
    % speed up computation in circulant case by smaller grid
    if min(abs(z_pts(jj)-z))<1.2
        RES2(jj) = svds(B - z_pts(jj)*speye(size(B)),1,'smallest');
    end
    parfor_progress(pf);
end

RES = reshape(RES,length(y_pts),length(x_pts));
RES2 = reshape(RES2,length(y_pts),length(x_pts));

%% Plot pseudospectra

v=(10.^(-6:0.5:0));

figure
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(max(real(RES),min(v))),log10(v),'LineColor','k',...
    'linewidth',0.5,'ShowText','off');
cbh=colorbar;
cbh.Ticks=log10(10.^(-8:1:0));
cbh.TickLabels=["1e-08","1e-07","1e-06","1e-05","1e-04","1e-03","1e-02","1e-01","1"];
clim([-5,0]);
colormap gray
hold on
plot(real(E),imag(E),'k.','markersize',12)
title('$\mathrm{Sp}_{\epsilon}(T_{50})$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
xlim([-12,17.01])
ylim([-16.01,17.01])
box on


figure
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(max(real(RES2),min(v))),log10(v),'LineColor','k',...
    'linewidth',0.5,'ShowText','off');
cbh=colorbar;
cbh.Ticks=log10(10.^(-8:1:0));
cbh.TickLabels=["1e-08","1e-07","1e-06","1e-05","1e-04","1e-03","1e-02","1e-01","1"];
clim([-5,0]);
colormap gray
hold on
plot(real(E2),imag(E2),'k.','markersize',12)
title('$\mathrm{Sp}_{\epsilon}(C_{50})$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
xlim([-12,17.01])
ylim([-16.01,17.01])
box on


