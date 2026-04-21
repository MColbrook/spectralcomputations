clear
close all

%% Construct matrix and implement finite section method

N = 50; % increase this to see convergence

a = ones(2*N+3,1);
b = 0*transpose((1).^(-(N+1):(N+1)));
H = spdiags([a,b,a],[-1,0,1],2*N+3,2*N+3);

nn = (-(N+1)):(N+1);
V = 1i + 0*nn;
V(abs(nn)>5)=0;

H2 = spdiags([a,transpose(V),a],[-1,0,1],2*N+3,2*N+3);
H0 = H2(2:end-1,2:end-1);
[V2,E] = eig(full(H0),'vector');

%% Compute pseudospectra using rectangular sections

xpts=-3:0.05:3;    ypts=-1:0.05:2; % increase resolution if desired
% RECOMMENDATION: when increasing resolution, do things adaptively after a
% coarse grid estimate of the spectrum
zpts=kron(xpts,ones(length(ypts),1))+1i*kron(ones(1,length(xpts)),ypts(:));    zpts=zpts(:);		% complex points where we compute pseudospectra
RES=0*zpts;
pf = parfor_progress(length(zpts));
pfcleanup = onCleanup(@() delete(pf));

II =speye(size(H2));

for jj=1:length(zpts)
    B = H2(:,2:end-1) - II(:,2:end-1)*zpts(jj);
    RES(jj) = svds(B,1,'smallest','MaxIterations',2000);
    parfor_progress(pf);
end
RES=reshape(RES,length(ypts),length(xpts));

%% Plot pseudospectra
v=(10.^(-2:0.25:0));
figure
hold on
contourf(reshape(real(zpts),length(ypts),length(xpts)),reshape(imag(zpts),length(ypts),length(xpts)),log10(max(abs(RES),min(v))),log10(v),'LineColor','k',...
    'linewidth',1,'ShowText','off');%,'FaceAlpha',0.6);
cbh=colorbar;
cbh.Ticks=log10(10.^(-4:1:0));
cbh.TickLabels=10.^(-4:1:0);
clim([-2,0])
colormap gray

hold on
plot(real(E),imag(E),'.r') % finite section eigenvalues


title('$\mathrm{Sp}_\epsilon(H)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)

ax=gca; ax.FontSize=18;
box on

%% CompSpec with wrong g
spec0 = cell(length(zpts),1);

pf = parfor_progress(length(zpts));
pfcleanup = onCleanup(@() delete(pf));
for jj=1:length(zpts)
    if RES(jj)<0.5
        I = find(abs(zpts-zpts(jj))<=RES(jj));
        s = zpts(I);
        d = RES(I);
        spec0{jj}=s(d==min(d));
    end
    parfor_progress(pf);
end
spec0=cell2mat(spec0);

%% CompSpec with correct g
spec = cell(length(zpts),1);
RES2 = RES*5;

pf = parfor_progress(length(zpts));
pfcleanup = onCleanup(@() delete(pf));
for jj=1:length(zpts)
    if RES2(jj)<0.5
        I = find(abs(zpts-zpts(jj))<=RES2(jj));
        s = zpts(I);
        d = RES2(I);
        spec{jj}=s(d==min(d));
    end
    parfor_progress(pf);
end

spec=cell2mat(spec);

%%

figure

tiledlayout(2,1,"TileSpacing","compact")
nexttile
spec0 = unique(round(spec0*100)/100);

plot(real(spec0),imag(spec0),'.k','markersize',10)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
title('Incorrect $g$','interpreter','latex','fontsize',18)
ax = gca; ax.FontSize = 18;

axis equal
axis([-2.5,2.5,-0.1,1.1])

nexttile

plot(real(spec),imag(spec),'.k','markersize',12)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
title('Correct $g$','interpreter','latex','fontsize',18)

axis equal
axis([-2.5,2.5,-0.1,1.1])
ax = gca; ax.FontSize = 18;



