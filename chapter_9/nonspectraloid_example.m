clear 
close all

%%%%%%%%%% CODE FOR NON-SPECTRALOID EXAMPLE (FIGURE 9.3) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compute pseudospectra and numerical range

B = [2,-1;1,-1];
E = eig(B);
p = NR(B,100); % compute numerical range

% compute pseudospectra
x_pts=-2:0.02:3;    y_pts=-1.5:0.02:1.5;
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);		% complex points where we compute pseudospectra
RES=0*z_pts;

pf = parfor_progress(length(z_pts));
pfcleanup = onCleanup(@() delete(pf));

for jj=1:length(z_pts)
    RES(jj) = svds(B - z_pts(jj)*eye(2),1,'smallest');
    parfor_progress(pf);
end
RES=reshape(RES,length(y_pts),length(x_pts));

%% Plot results

v = (10.^(-4:0.2:0));
figure
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(max(real(RES),min(v))),log10(v),'LineColor','k',...
    'linewidth',0.5,'ShowText','off');
cbh=colorbar;
cbh.Ticks=log10(10.^(-4:1:0));
cbh.TickLabels=10.^(-4:1:0);
clim([-4,0])
colormap gray
ax=gca; ax.FontSize=14;
axis tight equal;
hold on
plot(real(E),imag(E),'k.','markersize',5)
plot(real(p),imag(p),'-k','linewidth',2)
text(0,1.03,'$\partial W(B)$','interpreter','latex','fontsize',16,'Rotation',0,'background',[0.95,0.95,0.95])
title('$\mathrm{Sp}_{\epsilon}(B)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
annotation('doublearrow', [0.39 0.61], [0.5 0.5]-0.01,'linestyle',':','linewidth',2);
annotation('doublearrow', [0.39 0.66], [0.5 0.5]+0.03,'linestyle',':','linewidth',2);
text(0.4,-0.4,'$\rho(B)$','interpreter','latex','fontsize',16,'Rotation',0)
text(0.4,+0.32,'$\mu(B)$','interpreter','latex','fontsize',16,'Rotation',0)
ax=gca; ax.FontSize=18;
box on


function [p] = NR(A,k)
theta=(0:k)*2*pi/k;
p=0*theta;

for j=1:length(theta)
    H=exp(1i*theta(j))*A;
    H=(H+H')/2;
    [V,E]=eig(full(H));
    I=find(diag(E)==max(diag(E)));
    p(j)=(V(:,I(1)))'*(A*V(:,I(1)));
end
end