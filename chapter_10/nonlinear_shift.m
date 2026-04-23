clear
close all

%%%%%%%%%% CODE FOR NONLINEAR SHIFT EXAMPLE (FIGURE 10.2) %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up the parameters

f = @(z) sin(4*z).*(abs(z).^2+1);
N = 100;
S = spdiags(ones(2*N+3,1),1,2*N+3,2*N+3);% S = S';

xpts=-0.01:0.01:2;    ypts=-0.025:0.025:0.4;
zpts=kron(xpts,ones(length(ypts),1))+1i*kron(ones(1,length(xpts)),ypts(:));    zpts=zpts(:);		% complex points where we compute pseudospectra

RES=0*zpts+1;
RES2=RES;


%% Compute pseudospectra

% Sigma_1 algorithm
pf = parfor_progress(length(find(abs(abs(f(zpts))-1)<2)));
pfcleanup = onCleanup(@() delete(pf));
for jj=1:length(zpts)
    if abs(abs(f(zpts(jj)))-1)<2
        B = S-f(zpts(jj))*S';
        B = B(:,2:end-1);
        RES(jj) = svds(B,1,'smallest');
        parfor_progress(pf);
    end 
end

% finite section
pf = parfor_progress(length(zpts));
pfcleanup = onCleanup(@() delete(pf));
for jj=1:length(zpts)
    B = S-f(zpts(jj))*S';
    B = B(1:40,1:40); % finite section
    RES2(jj) = svds(B,1,'smallest');
    parfor_progress(pf);
end

RES=reshape(RES,length(ypts),length(xpts));
RES2=reshape(RES2,length(ypts),length(xpts));

%% Plot the results
close all
v=(10.^(-20:0.2:0));

figure % use fourfold symmetry
contourf(reshape(real(zpts),length(ypts),length(xpts)),reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES),min(v))),log10(v),'LineColor',[1,1,1]*0.3,...
    'linewidth',1,'linestyle','-','ShowText','off');
hold on
contourf(reshape(real(zpts),length(ypts),length(xpts)),-reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES),min(v))),log10(v),'LineColor',[1,1,1]*0.3,...
    'linewidth',1,'linestyle','-','ShowText','off');
contourf(-reshape(real(zpts),length(ypts),length(xpts)),reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES),min(v))),log10(v),'LineColor',[1,1,1]*0.3,...
    'linewidth',1,'linestyle','-','ShowText','off');
hold on
contourf(-reshape(real(zpts),length(ypts),length(xpts)),-reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES),min(v))),log10(v),'LineColor',[1,1,1]*0.3,...
    'linewidth',1,'linestyle','-','ShowText','off');
cbh=colorbar;
cbh.Ticks=log10(10.^(-20:1:0));
cbh.TickLabels=["1e-20","1e-19","1e-18","1e-17","1e-16","1e-15","1e-14","1e-13","1e-12","1e-11",...
    "1e-10","1e-9","1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1"];
clim([-2,0])
colormap gray
axis tight
title('Convergent Algorithm','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
box on
set(gca,'layer','top');


v=(10.^(-50:0.4*8:0));
figure
contourf(reshape(real(zpts),length(ypts),length(xpts)),reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES2),min(v))),log10(v),'LineColor',[1,1,1]*0.3,...
    'linewidth',1,'linestyle','-','ShowText','off');
hold on
contourf(reshape(real(zpts),length(ypts),length(xpts)),-reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES2),min(v))),log10(v),'LineColor',[1,1,1]*0.3,...
    'linewidth',1,'linestyle','-','ShowText','off');
contourf(-reshape(real(zpts),length(ypts),length(xpts)),reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES2),min(v))),log10(v),'LineColor',[1,1,1]*0.3,...
    'linewidth',1,'linestyle','-','ShowText','off');
hold on
contourf(-reshape(real(zpts),length(ypts),length(xpts)),-reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES2),min(v))),log10(v),'LineColor',[1,1,1]*0.3,...
    'linewidth',1,'linestyle','-','ShowText','off');
cbh=colorbar;
cbh.Ticks=log10(10.^(-30:5:0));
cbh.TickLabels=["1e-30","1e-25","1e-20","1e-15","1e-10","1e-5","1"];
clim([-30,0])
colormap gray
axis tight
title('Finite Section','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
box on
set(gca,'layer','top');


