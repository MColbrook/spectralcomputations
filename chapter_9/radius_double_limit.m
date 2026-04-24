clear
close all

%%%%%%%%%%%%%%%%%%%% CODE FOR EXAMPLE 9.1.12 %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot spectra

th = 0:0.01:2*pi*1.01;
z = exp(-1i*th);
a = z.^4-1i*z.^2;

figure
plot(real(a),imag(a),'k','linewidth',2)
title('$\mathrm{Sp}(A^2)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
axis equal
axis([-2,2,-2,2])

th = -pi/4:0.01:pi/4;
z = exp(-1i*th);
a = z.^4-1i*z.^2;

figure
a = sqrt(a);
plot(real(a),imag(a),'k','linewidth',2)
hold on
plot(-real(a),-imag(a),'k','linewidth',2)
plot(-real(a),imag(a),'k','linewidth',2)
plot(real(a),-imag(a),'k','linewidth',2)
title('$\mathrm{Sp}(A)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
axis equal
axis([-2,2,-2,2])


%% Pi_2 algorithm for spectral radius

xgrid = 1.4:0.001:1.44;
Nvec = round((10.^(1:0.05:4))/2);
egrid = zeros(length(Nvec),length(xgrid));

pf = parfor_progress(length(Nvec));
pfcleanup = onCleanup(@() delete(pf));
for ct = 1:length(Nvec)
    N = Nvec(ct)+2;
    a = repmat([1;1i;-1;-1i],N,1);
    b = 0*a + 1;
    A = spdiags([b,a],[2,1],4*N,4*N);
    B = speye(size(A));
    A = A(:,5:end-4);
    B = B(:,5:end-4);
    for j = 1:length(xgrid)
        C = A-xgrid(j)*1i*B; % use fact that maximum pseudospecral point occurs on imaginary axis
        egrid(ct,j) = svds(C,1,'smallest');
    end
    parfor_progress(pf);
end


n2 = [50,100,250,1000];
Gamma = zeros(length(Nvec),length(n2))+NaN;

for ii = 1: length(Nvec)
    for jj = 1:length(n2)
        I = find(egrid(ii,:)<1/n2(jj));
        if ~isempty(I)
            Gamma(ii,jj) = xgrid(max(I));
        end
    end
end

for jj = 1:4
    I1 = find(Gamma(:,jj)>1.4);
    Gamma(I1,jj)=smooth(Gamma(I1,jj));
end

figure
semilogx(Nvec*2,(Gamma(:,1)),'k','linewidth',3)
hold on
plot(Nvec*2,(Gamma(:,2)),':k','linewidth',2)
plot(Nvec*2,(Gamma(:,3)),'-.k','linewidth',2)
plot(Nvec*2,(Gamma(:,4)),'--k','linewidth',2)
plot(Nvec*2,Nvec*0+sqrt(2),'-k','linewidth',1)
ylim([1.4,1.44])
xlim([100,max(Nvec)*2])
xlabel('$n_1/2$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
legend({'$\epsilon=0.02$','$\epsilon=0.01$','$\epsilon=0.004$','$\epsilon=0.001$','$\rho(A)=\sqrt{2}$'},'interpreter','latex',...
    'fontsize',18,'location','southwest','NumColumns',2)







