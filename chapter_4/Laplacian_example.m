clear
close all
addpath(genpath('SpecSolve'));

%% Laplacian on periodic interval (done analytically)

f = chebfun(@(t) 1./(1.5+sin(t)),[-pi pi],'trig');
f = f/norm(f);
a = trigcoeffs(f)*sqrt(2*pi);


nn = 0:10;
I = (length(a)-1)/2+1;

J = [I-nn(2:end),I+nn];
y = a(J);
y = transpose(y(:));
x=[-nn(2:end),nn];

X1 = abs(x(:)).^2;
Y = abs(y(:)).^2;

x = [x;x;x];
y = [y;0*y+10^(-16);0*y+NaN];
x=x(:);
y=y(:);

y(x==0)=y(x==0)/sqrt(2);

figure
semilogy(abs(x).^2,2*abs(y).^2,'k','linewidth',2)
xlabel('$x$','interpreter','latex','fontsize',18)
title('Interval $[-\pi,\pi]_{\mathrm{per}}$','interpreter','latex','fontsize',18)
ylim([10^(-8),1])
xlim([-1,50])
ax=gca; ax.FontSize=18;
yticks(10.^(-8:2:0))

% smoothed measure

xx = -1:0.001:100;

epsilon = 0.1;
Y1 = sum(Y*epsilon^2/pi./(epsilon^2+(X1-xx).^2),1)*pi;
epsilon = 0.05;
Y2 = sum(Y*epsilon^2/pi./(epsilon^2+(X1-xx).^2),1)*pi;
epsilon = 0.01;
Y3 = sum(Y*epsilon^2/pi./(epsilon^2+(X1-xx).^2),1)*pi;

figure
semilogy(xx,Y1,':k','linewidth',2)
hold on
plot(xx,Y2,'--k','linewidth',2)
plot(xx,Y3,'k','linewidth',2)
xlabel('$x$','interpreter','latex','fontsize',18)
title('Interval $[-\pi,\pi]_{\mathrm{per}}$, $\pi\epsilon\cdot\mu_v^\epsilon(x)$','interpreter','latex','fontsize',18)
semilogy(abs(x).^2,2*abs(y).^2,'g','linewidth',1)
legend({'$\epsilon=0.1$','$\epsilon=0.05$','$\epsilon=0.01$'},'fontsize',16,'interpreter','latex','location','northeast')
ylim([10^(-8),1])
xlim([-1,50])
ax=gca; ax.FontSize=18;
yticks(10.^(-8:2:0))


%% Laplacian on real line (done computationally)

figure
xx = 0:0.001:1;
plot(xx,exp(-xx/2)./(2*sqrt(xx))/((sqrt(pi/2))),'k','linewidth',2)
ylim([0,5])
xlim([-0.1,0.501])
xlabel('$x$','interpreter','latex','fontsize',18)
title('Real line, $\rho_v(x)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
hold on
plot([0,0],[0,5],'--k')

X=-0.1:0.001:0.5;                                     %Evaluation pts
f=@(x) exp(-x.^2)/(sqrt(sqrt(pi/2)));               %Measure wrt f(r)
a={@(x) 0*x, @(x) 0*x, @(x) -1+0*x}; %Schrodinger op


mu1=diffMeas(a,f,X,0.1,'order',1);  %epsilon=0.1, m=1
mu2=diffMeas(a,f,X,0.05,'order',1);
mu3=diffMeas(a,f,X,0.01,'order',1);

figure
plot(X,mu1,':k','linewidth',2)
hold on
plot(X,mu2,'--k','linewidth',2)
plot(X,mu3,'k','linewidth',3)
plot(xx,exp(-xx/2)./(2*sqrt(xx))/((sqrt(pi/2))),'k','linewidth',1)
plot([0,0],[0,5],'--k')
ylim([0,5])
xlim([-0.1,0.501])

xlabel('$x$','interpreter','latex','fontsize',18)
title('Real line, $\mu_v^\epsilon(x)$','interpreter','latex','fontsize',18)
legend({'$\epsilon=0.1$','$\epsilon=0.05$','$\epsilon=0.01$','$\rho_v(x)$'},'fontsize',16,'interpreter','latex','location','northeast')
ax=gca; ax.FontSize=18;



