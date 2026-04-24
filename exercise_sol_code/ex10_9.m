% close all
clear

%% Set parameters
V = @(x) -5*exp(-abs(x)); % potential
N = 150; % 300
n = 200;
l = 500;


cntr = -(1+sqrt(3))/2;
Z0 = cntr;
LL = (sqrt(3)-1)/2.5;

contour = @(t) cntr + LL*cos(t)+LL*1i*sin(t);
tpts = linspace(0,2*pi,n+1)+2*pi/(2*n); tpts(end)=[];  % trap rule
quadpts = contour(tpts);
quadwts = (quadpts-cntr)/n;

%% Contour method
m = zeros(3,l);
N2 = [2,5,10];

for k = 1:length(N2)
    rng(1); wght = @(bb) abs(bb)<N2(k)+1;
    G = diag(wght(-N:N))*(randn(2*N+1,l)+1i*randn(2*N+1,l))/sqrt(2);
    
    for j = 1:n
        [U1,U2] = nlevp_solver(quadpts(j), G,V);
        U = quadwts(j)*U1 + conj(quadwts(j))*U2;
        m(k,:) = m(k,:) + dot(G,U)/2;
    end
end
m = real(m);

%%
figure
plot(1:l,cumsum(m(1,:).')'./(1:l),'k','linewidth',1)
hold on
plot(1:l,cumsum(m(2,:).')'./(1:l),'--k','linewidth',1)
plot(1:l,cumsum(m(3,:).')'./(1:l),':k','linewidth',1)
plot(1:l,0*(1:l)+2,'k','linewidth',2)
xlim([10,500])

xlabel('$K$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{tr}_K(S)$','interpreter','latex','fontsize',18)
legend({'$N_2=2$','$N_2=5$','$N_2=10$'},'location','east','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
box on
set(gca,'layer','top');
exportgraphics(gcf,'ex10_9_hutch.pdf','ContentType','vector','BackgroundColor','none','Colorspace','gray')




%% Code for solving linear systems

function [Sol1,Sol2] = nlevp_solver(z, RHS, V)
    N = (size(RHS,1)-1)/2;
    % Jacobi operator part
    S = spdiags(2-(-1).^((-(N+1):N+1)'),1,2*N+3,2*N+3); DIFF = (S+S')/2;
    
    % Potential part
    V1=sparse(diag(V(-(N+1):N+1)));
    V2=V1.^2;
    
    A0 = -DIFF+2*speye(size(DIFF))-V2;
    A0 = (A0+A0')/2;
    A1 = 2*V1;
    A2 = -speye(size(DIFF));

    bb = 1; % bandwidth

    A0 = A0(1+bb:end-bb,1+bb:end-bb);
    A1 = A1(1+bb:end-bb,1+bb:end-bb);
    A2 = A2(1+bb:end-bb,1+bb:end-bb);
    
    Sol1 = (A1+2*A2*z)*((A0+A1*z+A2*z^2)\RHS);
    Sol2 = ((A0+A1*z+A2*z^2)')\((A1+2*A2*z)'*RHS); % apply adjoint for symmetry trick


end


