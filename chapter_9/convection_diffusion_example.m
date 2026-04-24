clear
close all

%%%%%%% CODE FOR CONVECTION-DIFFUSION EXAMPLE (FIGURE 9.7) %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 10; % number of sine functions
nvec = [5,10,20,30];
tvec = 0:100;
pw = zeros(length(tvec),length(nvec));

for ii = 1:length(nvec)
    n=nvec(ii);
    F = chebfun(@(x) sin(pi/(2*n).*(1:N)*x)/sqrt(n),[0,2*n],'turbo'); % use Chebfun (very easy to code up!)
    L = chebop([0, 2*n]);
    L.op = @(u) diff(u, 2) + 2*diff(u) +0.5*u; 
    L.lbc = @(u) u;    L.rbc = @(u) u; % Dirichlet boundary conditions
    B=zeros(N,N,length(tvec));
    for jj=1:N
        u = expm(L, tvec, F(:,jj));       % exponential of the operator
        B(:,jj,:)=F'*u;
    end
    for jj=1:length(tvec)
        pw(jj,ii) = norm(squeeze(B(:,:,jj)));
    end
end


ff=figure;
semilogy(tvec,pw,'k','linewidth',2)
hold on

semilogy(tvec(tvec>60),10^20*exp(-tvec(tvec>60)*(0.5)),'--k','linewidth',1)
semilogy(tvec((tvec<30)&(tvec>0)),10*exp(tvec((tvec<30)&(tvec>0))*(0.5)),'--k','linewidth',1)

text(6,10000,'$\mathcal{O}(\exp(t/2))$','interpreter','latex','fontsize',16,'Rotation',30)
text(72,1000000,'$\mathcal{O}(\exp(-t/2))$','interpreter','latex','fontsize',16,'Rotation',-27)
text(35,1/1700000,'$n=5$','interpreter','latex','fontsize',16,'Rotation',-32)
text(48,1/40000,'$n=10$','interpreter','latex','fontsize',16,'Rotation',-29)
text(75,10,'$n=30$','interpreter','latex','fontsize',16,'Rotation',-27)
text(60,0.8,'$n=20$','interpreter','latex','fontsize',16,'Rotation',-27)
title('$\|\exp(tA_n)\|$','interpreter','latex','fontsize',18)

ylim([10^(-15),10^10])
xlabel('$t$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
ff.Position=[360.0000   97.6667  560.0000*1.2  420.0000];
