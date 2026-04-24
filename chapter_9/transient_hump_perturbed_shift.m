clear
close all

%%%%%%%%%%%%%% CODE FOR TRANSIENT HUMP (FIGURE 9.1) %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set the parameters

rng(0)
c = 1.01;

%% Build matrix

N = 2000;
S1 = [sparse(1,2*N+3); speye(2*N+2),sparse(2*N+2,1)]; S1 = S1';
A = S1;
A(N+2,N+2)=1; % rank-one perturbation
I = speye(2*N+3);

%% Pseudospec (this is used to compute pseudospectra abscissa)

z_pts=(c+0.01):0.001:1.2;
RES=0*z_pts;
pf = parfor_progress(length(z_pts));
pfcleanup = onCleanup(@() delete(pf));
for jj=1:length(z_pts)
    B1 = A - z_pts(jj)*I;
    B2 = B1';
    RES(jj) =  min(svds(B1(:,2:end-1),1,'smallest'),svds(B2(:,2:end-1),1,'smallest'));
    parfor_progress(pf);
end


%% Matrix powers (for larger powers it is more efficient to do in powers of two and combine)

A = A(2:end-1,2:end-1);

A = A/c;
B = A;

k = zeros(min(N-2,1000),1)+1;

for jj = 1:length(k)
    k(jj+1) = svds(B,1,'largest');
    B = B*A;
    plot((0:length(k)-1),k)
    pause(0.000001)
end

%% Plot the results

figure
plot((0:length(k)-1),k,'k','linewidth',2)
hold on
plot((0:length(k)-1),(0:length(k)-1)*0+max((z_pts-c)./RES),'--k','linewidth',1);
xlabel('$k$','interpreter','latex','fontsize',18)
title('$c^{-k}\|A^k\|$','interpreter','latex','fontsize',18)

ax=gca; ax.FontSize=18;



