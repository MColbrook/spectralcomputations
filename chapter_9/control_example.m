clear
close all

%%%%%%%%%%%%%%%%%% CODE FOR THE CONTROL EXAMPLE %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up the matrices

A = [5, -2;8 -4];
for j =1:12
    A = [A,sparse(size(A,1),size(A,2));sparse(size(A,1),size(A,2)),A];
end

N = size(A,1);
S1 = [sparse(1,2*N+3); speye(2*N+2),sparse(2*N+2,1)]; S1 = S1'+S1;
A = -S1(1:N,1:N)*A-speye(N,N);
N = 300; % number of columns used
Aadj = A';


%% Compute pseudospectra of finite sections

xpts = -2:0.05:0.5;    ypts = -0.01:0.05:1.5; % increse resolution if wanted
zpts=kron(xpts,ones(length(ypts),1))+1i*kron(ones(1,length(xpts)),ypts(:));    zpts=zpts(:);		% complex points where we compute pseudospectra
RES=0*zpts;

pf = parfor_progress(length(zpts));
pfcleanup = onCleanup(@() delete(pf));
V = rand(N,1); U = rand(N,1);
for jj=1:length(zpts)
    B = A(1:N,1:N) - speye(N,N)*zpts(jj);
    RES(jj)= svds(B,1,'smallest','MaxIterations',5000);%,'RightStartVector',V);%,'LeftStartVector',U);
    parfor_progress(pf);
end

RES=reshape(RES,length(ypts),length(xpts));


%% Compute pseudospectra using rectangular sections

RES2=0*zpts;

pf = parfor_progress(length(zpts));
pfcleanup = onCleanup(@() delete(pf));
for jj=1:length(zpts)
    B = A(1:(N+4),1:N) - speye(N+4,N)*zpts(jj);
    RES2(jj) = svds(B,1,'smallest','MaxIterations',5000);
    B = Aadj(1:(N+4),1:N) - speye(N+4,N)*conj(zpts(jj));
    RES2(jj) = min(RES2(jj),svds(B,1,'smallest','MaxIterations',5000));
    parfor_progress(pf);
end
RES2=reshape(RES2,length(ypts),length(xpts));
RES2(isnan(RES2))=10;


%% Plot results
close all

E = eig(full(A(1:N,1:N)));
v=(10.^(-16:0.33333:0));
figure
contourf(reshape(real(zpts),length(ypts),length(xpts)),reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES),min(v))),log10(v),'LineColor','k',...
    'linewidth',1,'ShowText','off');
hold on
contourf(reshape(real(zpts),length(ypts),length(xpts)),-reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES),min(v))),log10(v),'LineColor','k',...
    'linewidth',1,'ShowText','off');
cbh=colorbar;
cbh.Ticks=log10(10.^(-16:2:0));
cbh.TickLabels=["1e-16","1e-14","1e-12","1e-10","1e-08","1e-06","1e-04","1e-02","1"];
clim([-9,0])
colormap gray
ax=gca; ax.FontSize=14;
axis tight

axis([min(xpts),max(xpts),-max(ypts),max(ypts)])
title('Finite Section','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
box on


figure
contourf(reshape(real(zpts),length(ypts),length(xpts)),reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES2),min(v))),log10(v),'LineColor','k',...
    'linewidth',1,'ShowText','off');
hold on
contourf(reshape(real(zpts),length(ypts),length(xpts)),-reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES2),min(v))),log10(v),'LineColor','k',...
    'linewidth',1,'ShowText','off');
cbh=colorbar;
cbh.Ticks=log10(10.^(-16:2:0));
cbh.TickLabels=["1e-16","1e-14","1e-12","1e-10","1e-08","1e-06","1e-04","1e-02","1"];
clim([-9,0])
colormap gray
ax=gca; ax.FontSize=14;
axis tight
hold on
axis([min(xpts),max(xpts),-max(ypts),max(ypts)])
title('\texttt{PseudoSpec}','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)

ax=gca; ax.FontSize=18;
box on
