clear
close all

%%%%%%%%%%%%%%%% CODE FOR NH ANDERSON MODEL EXAMPLE %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set paramaters
n = 500; % truncation size is 2n+1. Increase n for convergence. Use parallelisation for very large n.
tau = 0.5; % coupling parameter
p = 0.5; % Bernoulli parameter

x_pts = -4:0.05:0.05;
y_pts = -1.4:0.025:0.025;

%% Construct matrix

S1 = [sparse(1,2*n+3); speye(2*n+2),sparse(2*n+2,1)];
S2 = S1';
H = exp(tau)*S1 + exp(-tau)*S2 + spdiags(sign((rand(2*n+3,1)-p)),0,2*n+3,2*n+3);
I = speye(2*n+3);

%% Pseudospec

z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);		% complex points where we compute pseudospectra
RES=0*z_pts;
RES_fs=0*z_pts;
RES_per=0*z_pts;

pf = parfor_progress(length(z_pts));
pfcleanup = onCleanup(@() delete(pf));

for jj=1:length(z_pts)
    B1 = H - z_pts(jj)*I;
    B2 = B1';
    RES(jj) =  min(svds(B1(:,2:end-1),1,'smallest'),svds(B2(:,2:end-1),1,'smallest'));
    B1 = B1(2:end-1,2:end-1);
    B2 = B2(2:end-1,2:end-1);
    RES_fs(jj) =  svds(B1,1,'smallest');
    B1(1,end) = exp(tau); % periodic BCs
    B1(end,1) = exp(-tau);
    RES_per(jj) =  svds(B1,1,'smallest');
    parfor_progress(pf);
end
RES=reshape(RES,length(y_pts),length(x_pts));
RES_fs=reshape(RES_fs,length(y_pts),length(x_pts));
RES_per=reshape(RES_per,length(y_pts),length(x_pts));

%%
v=(10.^(-8:0.5:0));

figure
contourf(x_pts,y_pts,log10(max(real(RES),min(v))),log10(v),'LineColor','k','linewidth',1,'ShowText','off');
hold on % use symmetry of large n limit
contourf(-x_pts,y_pts,log10(max(real(RES),min(v))),log10(v),'LineColor','k','linewidth',1,'ShowText','off');
contourf(x_pts,-y_pts,log10(max(real(RES),min(v))),log10(v),'LineColor','k','linewidth',1,'ShowText','off');
contourf(-x_pts,-y_pts,log10(max(real(RES),min(v))),log10(v),'LineColor','k','linewidth',1,'ShowText','off');

cbh=colorbar;
cbh.Ticks=log10(10.^(-8:1:0));
cbh.TickLabels=["1e-08","1e-07","1e-06","1e-05","1e-04","1e-03","1e-02","1e-01","1"];
clim([-6,0])

colormap gray
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',13)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',13)
title('Output of \texttt{PseudoSpec} ($p=1/2$)','interpreter','latex','fontsize',13)
ax=gca; ax.FontSize=14;
axis tight equal;
set(gca,'ydir','normal');
axis([-4,4,-1.4,1.4])
box on


figure
contourf(x_pts,y_pts,log10(max(real(RES_fs),min(v))),log10(v),'LineColor','k','linewidth',1,'ShowText','off');
hold on % use symmetry
contourf(-x_pts,y_pts,log10(max(real(RES_fs),min(v))),log10(v),'LineColor','k','linewidth',1,'ShowText','off');
contourf(x_pts,-y_pts,log10(max(real(RES_fs),min(v))),log10(v),'LineColor','k','linewidth',1,'ShowText','off');
contourf(-x_pts,-y_pts,log10(max(real(RES_fs),min(v))),log10(v),'LineColor','k','linewidth',1,'ShowText','off');

cbh=colorbar;
cbh.Ticks=log10(10.^(-8:1:0));
cbh.TickLabels=["1e-08","1e-07","1e-06","1e-05","1e-04","1e-03","1e-02","1e-01","1"];
clim([-6,0])

colormap gray
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',13)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',13)
title('Finite Section with no B.C.s ($p=1/2$)','interpreter','latex','fontsize',13)
ax=gca; ax.FontSize=14;
axis tight equal;
set(gca,'ydir','normal');
axis([-4,4,-1.4,1.4])
box on

figure
contourf(x_pts,y_pts,log10(max(real(RES_per),min(v))),log10(v),'LineColor','k','linewidth',1,'ShowText','off');
hold on % use symmetry
contourf(-x_pts,y_pts,log10(max(real(RES_per),min(v))),log10(v),'LineColor','k','linewidth',1,'ShowText','off');
contourf(x_pts,-y_pts,log10(max(real(RES_per),min(v))),log10(v),'LineColor','k','linewidth',1,'ShowText','off');
contourf(-x_pts,-y_pts,log10(max(real(RES_per),min(v))),log10(v),'LineColor','k','linewidth',1,'ShowText','off');

cbh=colorbar;
cbh.Ticks=log10(10.^(-8:1:0));
cbh.TickLabels=["1e-08","1e-07","1e-06","1e-05","1e-04","1e-03","1e-02","1e-01","1"];
clim([-6,0])

colormap gray
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',13)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',13)
title('Finite Section with periodic B.C.s ($p=1/2$)','interpreter','latex','fontsize',13)
ax=gca; ax.FontSize=14;
axis tight equal;
set(gca,'ydir','normal');
axis([-4,4,-1.4,1.4])
box on

