close all
clear

%%%%%%%%%%% CODE FOR KLEIN-GORDON EXAMPLE (FIGURE 10.3) %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set parameters
V = @(x) -5*exp(-abs(x)); % potential
N = 100; 

%% Build operator

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


%% Finite section poly eigenvalue problem

[X,e] = polyeig(A0(1+bb:end-bb,1+bb:end-bb),A1(1+bb:end-bb,1+bb:end-bb),A2(1+bb:end-bb,1+bb:end-bb));
R =  A0(:,1+bb:end-bb)*X + A1(:,1+bb:end-bb)*X.*transpose(e) + A2(:,1+bb:end-bb)*X.*transpose(e.^2);
R = vecnorm(R)./vecnorm(X); % residual


%% Pseudospectra using convergent algorithm

xpts=-3.7:0.02:2.1;    ypts=-.02:0.02:0.8;
zpts=kron(xpts,ones(length(ypts),1))+1i*kron(ones(1,length(xpts)),ypts(:));    zpts=zpts(:);		% complex points where we compute pseudospectra

RES=0*zpts+1;

pf = parfor_progress(length(zpts));
pfcleanup = onCleanup(@() delete(pf));

for jj=1:length(zpts)
    B = A0 + A1*zpts(jj) + A2*zpts(jj)^2;
    B = B(:,1+bb:end-bb);
    RES(jj) = svds(B,1,'smallest','Tolerance',1e-5,'MaxIterations',10000);
    parfor_progress(pf);  
end

RES=reshape(RES,length(ypts),length(xpts));

%% Plot the results


close all
v=(10.^(-20:0.2:2));
f=figure;
contourf(reshape(real(zpts),length(ypts),length(xpts)),reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES),min(v))),log10(v),'LineColor',[1,1,1]*0.3,...
    'linewidth',1,'linestyle','-','ShowText','off');
hold on
contourf(reshape(real(zpts),length(ypts),length(xpts)),-reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES),min(v))),log10(v),'LineColor',[1,1,1]*0.3,...
    'linewidth',1,'linestyle','-','ShowText','off');
cbh=colorbar;
cbh.Ticks=log10(10.^(-20:1:0));
cbh.TickLabels=["1e-20","1e-19","1e-18","1e-17","1e-16","1e-15","1e-14","1e-13","1e-12","1e-11",...
    "1e-10","1e-9","1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1"];
clim([-3,0])
colormap gray
axis tight
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
box on
set(gca,'layer','top');
hold on
plot(real(e),imag(e),'.r','markersize',16)

set(gca,'layer','top');
f.Position=[160.0000   97.6667  560.0000*2  420.0000];
axis([min(xpts(:)),max(xpts(:)),-max(ypts(:)),max(ypts(:))])

x1 = [-1.415,-1.415];
y1 = -[0.9 0.05];

quiver( x1(1),y1(1),x1(2)-x1(1),y1(2)-y1(1),0,'m','linewidth',4,'MaxHeadSize' ,0.8) % arrow showing pollution

x1 = x1+2.82;
y1 = -[0.9 0.05];
quiver( x1(1),y1(1),x1(2)-x1(1),y1(2)-y1(1),0,'m','linewidth',4,'MaxHeadSize' ,0.8) % arrow showing pollution





