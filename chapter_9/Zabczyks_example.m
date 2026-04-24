close all
clear

%%%%%%%%%%%%% CODE FOR ZABCZYKS EXAMPLE (FIGURE 9.8) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J1 = 0; J1 = J1+1i*eye(1);
J2 = eye(3); J2(:,end)=[]; J2(1,:)=[]; J2 = J2+2i*eye(2);
J3 = eye(4); J3(:,end)=[]; J3(1,:)=[]; J3 = J3+3i*eye(3);
J4 = eye(5); J4(:,end)=[]; J4(1,:)=[]; J4 = J4+4i*eye(4);
J5 = eye(6); J5(:,end)=[]; J5(1,:)=[]; J5 = J5+5i*eye(5);
J6 = eye(7); J6(:,end)=[]; J6(1,:)=[]; J6 = J6+6i*eye(6);
J7 = eye(8); J7(:,end)=[]; J7(1,:)=[]; J7 = J7+7i*eye(7);
J8 = eye(9); J8(:,end)=[]; J8(1,:)=[]; J8 = J8+8i*eye(8);
J9 = eye(10); J9(:,end)=[]; J9(1,:)=[]; J9 = J9+9i*eye(9);
A = sparse(blkdiag(J1,J2,J3,J4,J5,J6,J7,J8,J9));

x_pts=-0.03:0.03:2;    y_pts=0:0.03:9; % use symmetry to reduce grid size
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);		% complex points where we compute pseudospectra
RES=0*z_pts;

pf = parfor_progress(length(z_pts));
pfcleanup = onCleanup(@() delete(pf));

for jj=1:length(z_pts)
    B = A - z_pts(jj)*speye(size(A));
    RES(jj) = svds(B,1,'smallest');
    parfor_progress(pf);
end
RES=reshape(RES,length(y_pts),length(x_pts));


v=(10.^(-4:0.4:0));
figure
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(max(real(RES),min(v))),log10(v),'LineColor','k',...
    'linewidth',0.5,'ShowText','off');
contourf(reshape(-real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(max(real(RES),min(v))),log10(v),'LineColor','k',...
    'linewidth',0.5,'ShowText','off'); % use symmetry to reduce grid size
cbh=colorbar;
cbh.Ticks=log10(10.^(-4:1:0));
cbh.TickLabels=10.^(-4:1:0);
clim([-4,0])
colormap gray
axis equal;
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
box on




