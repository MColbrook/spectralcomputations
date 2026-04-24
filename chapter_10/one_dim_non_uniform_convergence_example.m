clear
close all

%%%%%%%%%%%%%%%%%%%%%% CODE FOR FIGURE 10.1 %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon = 0.25;

xpts=-1:0.002:1;    ypts=-1:0.002:1;
zpts=kron(xpts,ones(length(ypts),1))+1i*kron(ones(1,length(xpts)),ypts(:));    zpts=zpts(:);
R1 = 0.00+abs(zpts.*(1-zpts))<=epsilon; % test for inclusion in pseudospectra
R2 = 0.00+abs(zpts)<1.1;
R3 = R2; R3(abs(zpts)>=1)=20; % to exclude points outside the disc

RES=reshape(1-R1.*R2.*R3,length(ypts),length(xpts));

v = (10.^(-1:0));
cSCALE = [-2,-0.5]; % scale for logarithmic epsilon

figure
contourf(reshape(real(zpts),length(ypts),length(xpts)),reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES),min(v))),log10(v),'LineColor',[1,1,1]*0,...
    'linewidth',1.5,'linestyle','-','ShowText','off');
clim(cSCALE)
colormap gray
set(gca,'Color',[1,1,1])
ax=gca; ax.FontSize=14;
hold on
r = 1;                % radius of disk
theta = linspace(0,2*pi,500);
x = r*cos(theta);
y = r*sin(theta);

% Outer rectangle covering the plotting region
Xrect = [-1  1  1 -1];  
Yrect = [-1 -1  1  1];

% Combine rectangle and disk to make a polygon with a hole
pgon = polyshape({Xrect, x}, {Yrect, y}, 'Simplify', false);
plot(pgon, 'FaceColor', [1 1 1], 'EdgeColor', 'none','FaceAlpha',0); 

r = 1;                % radius of disk
theta = linspace(0,2*pi,500);
x = r*cos(theta);
y = r*sin(theta);

% Outer rectangle covering the plotting region
Xrect = [-1  1  1 -1]*1.1;  
Yrect = [-1 -1  1  1]*1.1;

% Combine rectangle and disk to make a polygon with a hole
pgon = polyshape({Xrect, x}, {Yrect, y}, 'Simplify', false);
plot(pgon, 'FaceColor', [1 1 1], 'EdgeColor', 'none','FaceAlpha',1); 


plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'--k','linewidth',1.5)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
title(sprintf('$\\mathrm{Sp}_{\\epsilon}(T)$, $\\epsilon=%.2f$',epsilon),'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

axis equal
box on
ylim([-1,1])
xlim([-1,1])




