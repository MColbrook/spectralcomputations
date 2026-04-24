clear
close all

% Add utils and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'FEM_data');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%% CODE FOR MAXWELL EIG PROBLEM (2D) %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Spectral pollution without operator folding
N = 8;
order = 1;
T = sprintf('square_FEM_N%d_order%d.mat',N,order); % FEM matrices computed using DOLFIN
load(T,'B1','A1');

% enforce boundary condition
v = vecnorm(B1);
I = find(v>10^(-12));
B1 = B1(I,I);
A1 = A1(I,I);

E=eigs(A1,B1,100,-1);
E=sort(E);

f = figure;
plot(1:min(50,length(E)),E(1:min(50,length(E))),'k.','markersize',20);
ylim([0,10])
xlim([0,50])
hold on
plot([0,50],[0,0],'--k','linewidth',1)
plot([0,50],[1,1],'--k','linewidth',1)
plot([0,50],[2,2],'--k','linewidth',1)
plot([0,50],[4,4],'--k','linewidth',1)
plot([0,50],[5,5],'--k','linewidth',1)
plot([0,50],[8,8],'--k','linewidth',1)
plot([0,50],[9,9],'--k','linewidth',1)
plot([0,50],[10,10],'--k','linewidth',1)
ax=gca; ax.FontSize=18;
f.Position = 1.0e+02 *[3.600000000000000   0.976666666666667   4.203333333333332   4.200000000000000];
title(sprintf('$N=%d$',N),'interpreter','latex','fontsize',18)
xlabel('Eigenvalue Number','interpreter','latex','fontsize',18)
ylabel('$\omega^2$','interpreter','latex','fontsize',18)
box off


%% Distspec plot
clear

N = 32;
order = 1;

% import FEM matrices
T = sprintf('square_FEM_N%d_order%d.mat',N,order);
load(T);

J = 1:size(A,1);
J2 = [LD1(:);LD2(:);LD3(:);LD4(:)]+1; % +1 due to python starting at 0
J(J2) = [];

I = J;
A = A(I,I); A = (A+A')/2;
A0 = 1i*A0(I,I); A0 = (A0+A0')/2;
B = B(I,I); B = (B+B')/2;

% enforce boundary condition
v = vecnorm(B1);
I = find(v>10^(-12));
B1 = B1(I,I);
A1 = A1(I,I);

% compute approximation of distance to the spectrum
X = 0.5:0.01:4;
dist = 0*X;

pf = parfor_progress(length(X));
pfcleanup = onCleanup(@() delete(pf));
for jj=1:length(X)
    M1 = A-X(jj)*(A0'+A0)+abs(X(jj))^2*B;
    M2 = B;
    dist(jj) = min(sqrt(eigs(M1,M2,1,'smallestabs','maxit',1000)));
    parfor_progress(pf);
end


figure
plot(X,dist,'k','linewidth',1)
hold on

nn = [0 1 2 4 5 8 9 10];
plot(sqrt(nn),0*nn,'.k','markersize',20)
ax=gca; ax.FontSize=18;
xlabel('$x$','interpreter','latex','fontsize',18)
title(sprintf('$N=%d$',N),'interpreter','latex','fontsize',18)

legend({'$\Phi_{N}(x,M_{\mathrm{2D}})$','Exact $\omega$'},'fontsize',16,'interpreter','latex','location','northeast')

xlim([0.5,sqrt(10)])


%% Convergence plot
clear
EE = [1,2,4,5,8,9];
mult = [2,1,2,2,1,2];

Nvec = 2.^(1:7);

Ecomp = zeros(sum(mult),length(Nvec));
Err = zeros(sum(mult),length(Nvec));
ErrB = zeros(sum(mult),length(Nvec));
DOF = zeros(length(Nvec),1);

ct = 1;
for N = Nvec
    T = sprintf('square_FEM_N%d_order2.mat',N);
    load(T);

    J = 1:size(A,1);
    J2 = [LD1(:);LD2(:);LD3(:);LD4(:)]+1; % +1 due to python starting at 0
    J(J2) = [];
    
    I = J;
    A = A(I,I); A = (A+A')/2;
    A0 = 1i*A0(I,I); A0 = (A0+A0')/2;
    B = B(I,I); B = (B+B')/2;
    
    % enforce boundary condition
    v = vecnorm(B1);
    I = find(v>10^(-12));
    B1 = B1(I,I);
    A1 = A1(I,I);

    DOF(ct) = size(A,1);
    
    for jj = 1:length(EE)
        x = sqrt(EE(jj)+0.1); % in practice take a point near enough to the eigenvalue
        M1 = A-x*(A0'+A0)+abs(x)^2*B;
        M2 = B;
        [V,~] = eigs(M1,M2,mult(jj),0);
        RHO = real(dot(V,A0*V));
        EP = real(dot(V,A*V-(A0'+A0)*V.*RHO+B*V.*RHO.^2));

        EEE = EE;
        EEE(jj) = [];
        EEE = sqrt([0;EEE(:);10]);
        delta = 0*RHO;
        for jjj=1:length(RHO)
            delta(jjj) = min(abs(RHO(jjj)-EEE)); % in practice, compute lower bound on delta using distspec
        end

        
        EP = EP./delta;

        if jj>1
            II = (sum(mult(1:jj-1))+1):sum(mult(1:jj));
        else
            II = 1:sum(mult(1:jj));
        end
        Ecomp(II,ct) = RHO.^2;
        Err(II,ct) = EP;
        ErrB(II,ct) = abs(RHO-sqrt(EE(jj)));
    end
    ct = ct+1;
end


%% Plot convergence results
Err(Err>10)=NaN;
figure
loglog(Nvec,Err,'.-k','linewidth',1,'markersize',15)
xlim([1,128*2])
ax=gca; ax.FontSize=18;

hold on
% plot(Nvec(2:end-3),1./Nvec(2:end-3).^2/10,'--k','linewidth',2)
plot(Nvec(2:end-3),1./Nvec(2:end-3).^4/1000,'--k','linewidth',2)
grid minor
xlabel('$N$','interpreter','latex','fontsize',18)
% return

% text(4.1,10^(-3)*1.5,'$\mathcal{O}(N^{-2})$','interpreter','latex','fontsize',18,'Rotation',-34)
% title('A Posteriori Bound (Lagrange Order 1)','interpreter','latex','fontsize',18)
% ylim([0.0001,10])


text(4.1,10^(-7)*2,'$\mathcal{O}(N^{-4})$','interpreter','latex','fontsize',18,'Rotation',-32)
title('A Posteriori Bound (Lagrange Order 2)','interpreter','latex','fontsize',18)





