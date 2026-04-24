close all
clear

%%%%% CODE FOR KLEIN-GORDON CONTOUR EXAMPLE (SECTION 10.3.3) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set parameters

rng(0)
V = @(x) -5*exp(-abs(x)); % potential
N = 300; 
n = 200;
l = 8;
r = 8;
p = 1;


cntr = -(1+sqrt(3))/2;
Z0 = cntr;
LL = (sqrt(3)-1)/2.5;

contour = @(t) cntr + LL*cos(t)+LL*1i*sin(t);
jacobian = @(t) -LL*sin(t)+LL*1i*cos(t); 
tpts = linspace(0,2*pi,n+1)+2*pi/(2*n); tpts(end)=[];  % trap rule
jacobianpts=jacobian(tpts);
quadpts = contour(tpts);
quadwts = (2*pi)/n * ones(1,n);



wght = @(bb) 0*bb +1; % no decay (does not avoid spectral pollution!)
TT=sprintf('No decay for $\\{c_\\beta\\}$');

% wght = @(bb) exp(-abs(bb).^2/1000); % decying weight (works!)
% TT=sprintf('Decay for $\\{c_\\beta\\}$');



%% Contour method

% random probing matrices
rng(0);
F = diag(wght(-N:N))*(randn(2*N+1,l)+1i*randn(2*N+1,l))/sqrt(2);
G = diag(wght(-N:N))*(randn(2*N+1,l)+1i*randn(2*N+1,l))/sqrt(2);


A = cell({2*p});

for jj = 0:2*p-1
    A{jj+1} = zeros(2*N+1,r);
end

for j = 1:n
    U = nlevp_solver(quadpts(j), G,V);
    for jj = 0:2*p-1
        A{jj+1} = A{jj+1} + quadwts(j)*jacobianpts(j)*(quadpts(j)-Z0)^jj*U/2i/pi;
    end
end

B0 = zeros(p*l,p*r); B1 = B0;
for aa = 1:p
    for bb = 1:p
        B0((l*(aa-1)+1):aa*l,(r*(bb-1)+1):bb*r) = F'*A{aa+bb-1};
        B1((l*(aa-1)+1):aa*l,(r*(bb-1)+1):bb*r) = F'*A{aa+bb};
    end
end

C = [];
for bb = 1:p
    C = [C,A{bb}];
end


% Reduced SVD of A0:
[V0, S0, W0] = svd(B0, 0);


% Figure 10.9
figure
semilogy(diag(S0)/S0(1,1),'k.--','linewidth',2,'markersize',30)
grid on
hold on
ylim([10^(-18),1])
xlabel('$j$','interpreter','latex','fontsize',18)
title(TT,'interpreter','latex','fontsize',18)
ylabel('$\sigma_j(\tilde{B}_0)/\sigma_1(\tilde{B}_0)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=17;


%% Plot spurious eigenvector (when c_n=1)

[VV, lam] = eig(V0(:,1:3)'*B1*W0(:,1:3),S0(1:3,1:3),'vector');
[~,I] = sort(real(lam),'ascend'); VV = VV(:,I); lam = lam(I);
V = C*W0(:,1:3)*VV; V = V./vecnorm(V);

figure
semilogy(-N:N,abs(V(:,2)),'.k','markersize',16)
grid on
ax=gca; ax.FontSize=18;
xlabel('$\beta$','interpreter','latex','fontsize',24)
ylabel('$|d_{\beta}^{(N)}|$','interpreter','latex','fontsize',24)
title('$N=300$','interpreter','latex','fontsize',24)
xlim([-N,N])




%% Code for solving linear systems

function Sol = nlevp_solver(z, RHS, V)
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
    
    Sol = (A0+A1*z+A2*z^2)\RHS;
end


