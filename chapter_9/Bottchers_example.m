clear
close all

%%%%%%%%%%%% CODE FOR BOTTCHER'S EXAMPLE (FIGURE 9.2) %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Spectrum plot

A = spdiags(ones(50,2),[-1,2],50,50)*10/19;

figure
E=eig(full(A));
plot(real(E),imag(E),'.k','markersize',14)
th=0:0.001:2*pi;
a=(exp(-1i*th)+exp(2i*th))*10/19;
hold on
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
plot(real(a),imag(a),'-k','linewidth',1)
hold on
plot(cos(th)*15*2^(1/3)/19,sin(th)*15*2^(1/3)/19,'--k','linewidth',1)
ax=gca; ax.FontSize=18;
legend({'$\mathrm{Sp}(A_{50})$','$a(\bf{T})$','$|z|=\frac{15}{19}\sqrt[3]{2}$'},'interpreter','latex','fontsize',16,'location','southeast')
axis  equal
box on
axis([-1.1,1.1,-1.1,1.1])


%% Noncommuting limits

A1 = spdiags(ones(100,2),[-1,2],100,100)*10/19; B1 = A1;
A2 = spdiags(ones(200,2),[-1,2],200,200)*10/19; B2 = A2;
A3 = spdiags(ones(300,2),[-1,2],300,300)*10/19; B3 = A3;
A4 = spdiags(ones(400,2),[-1,2],400,400)*10/19; B4 = A4;

pw = zeros(5000,4);

pf = parfor_progress(5000);
pfcleanup = onCleanup(@() delete(pf));
for k = 1:5000 % one can speed this up by computing A^{2^j} for different j, but this is sufficiently fast for this example
    pw(k,1) = svds(A1,1,'largest'); A1 = A1*B1;
    pw(k,2) = svds(A2,1,'largest'); A2 = A2*B2;
    pw(k,3) = svds(A3,1,'largest'); A3 = A3*B3;
    pw(k,4) = svds(A4,1,'largest'); A4 = A4*B4;
    parfor_progress(pf);
end

figure;
semilogy(1:5000,pw(:,1:4),'k','linewidth',2)
hold on
tvec=0:5000;
semilogy(tvec((tvec>1000)&(tvec<3000)),(15*2^(1/3)/19).^tvec((tvec>1000)&(tvec<3000)),'--k','linewidth',1)
semilogy(tvec((tvec<2000)&(tvec>0)),1000*(20/19).^tvec((tvec<2000)&(tvec>0)),'--k','linewidth',1)
ylim([10^(-10),10^30])
text(300,10^14,'$\mathcal{O}((20/19)^k)$','interpreter','latex','fontsize',16,'Rotation',64)
text(1400,10^(-1),'$\mathcal{O}((15\sqrt[3]{2}/19)^k)$','interpreter','latex','fontsize',16,'Rotation',-10)
text(4000,10^(-2),'$n=100$','interpreter','latex','fontsize',16,'Rotation',-10)
text(4000,10^(8),'$n=200$','interpreter','latex','fontsize',16,'Rotation',-8)
text(4000,10^(17),'$n=300$','interpreter','latex','fontsize',16,'Rotation',-8)
text(4000,10^(25.5),'$n=400$','interpreter','latex','fontsize',16,'Rotation',-8)
title('$\|A_n^k\|$','interpreter','latex','fontsize',18)
xlabel('$k$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;



