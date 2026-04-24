clear
close all



%% Compute essential pseudospectra
N = 1000;
N2 = 2;
xpts=-1:0.025:1;    ypts=-1:0.025:0.25;
zpts=kron(xpts,ones(length(ypts),1))+1i*kron(ones(1,length(xpts)),ypts(:));    zpts=zpts(:);		% complex points where we compute pseudospectra

RES1=0*zpts+1; RES2=RES1;
pf = parfor_progress(length(zpts));
pfcleanup = onCleanup(@() delete(pf));

for jj=1:length(zpts)
    L = mat_setup(zpts(jj),N,[]);
    RES1(jj) = min(svd(L(:,N2:end)));
    L = mat_setup2(zpts(jj),N,[]);
    RES2(jj) = min(svd(L(:,N2:end)));
    parfor_progress(pf);  
end



%% Plot the results
RES=reshape(min(RES1,RES2),length(ypts),length(xpts));

v=(10.^(-3:0.15:0));
figure
contourf(reshape(real(zpts),length(ypts),length(xpts)),reshape(imag(zpts),length(ypts),length(xpts)),log10(max(real(RES),min(v))),log10(v),'LineColor',[1,1,1]*0,...
    'linewidth',1,'linestyle','-','ShowText','off');
cbh=colorbar;
cbh.Ticks=log10(10.^(-20:1:0));
cbh.TickLabels=["1e-20","1e-19","1e-18","1e-17","1e-16","1e-15","1e-14","1e-13","1e-12","1e-11",...
    "1e-10","1e-9","1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1"];
clim([-3,0])
colormap gray
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
box on
set(gca,'layer','top');

exportgraphics(gcf,'ex10_8_pseudospectra.pdf','ContentType','vector','BackgroundColor','none','Colorspace','gray')






function [L,Z] = mat_setup(lam,N,X)

%% Set up matrices using Laguerre functions

NN = N+1;
D = zeros(NN,NN);
for jj = 1:NN
    D(1:jj-1,jj)= -1;
    D(jj,jj) = -1/2;
end

L = D*D + lam^2*eye(NN);

% Basis recombination

B = [zeros(1,NN-1);eye(NN-1,NN-1)];
for jj = 1:NN-1
     B(1,jj) = (-(jj+1/2)-1i*lam)/(0.5+1i*lam);
end


[Q,~] = qr(B,"econ");

L = L*Q;
NN = NN - 1;
L = L(:,1:NN);


if ~isempty(X)
    X=X(:);
    LL = zeros(length(X),N+1);
    LL(:,1) = exp(-X/2);
    LL(:,2) = (1 - X).*exp(-X/2);
    for jj=2:N
        k = jj - 1;
        LL(:,jj+1) = (2*k+1 - X)./(k+1).*LL(:,jj)-k/(k+1)*LL(:,jj-1);
    end
    Z = LL*Q;
else
    Z = [];
end

end




function [L,Z] = mat_setup2(lam,N,X)

%% Set up matrices using Laguerre functions

NN = N+1;
D = zeros(NN,NN);
for jj = 1:NN
    D(1:jj-1,jj)= -1;
    D(jj,jj) = -1/2;
end

L = D*D + conj(lam)^2*eye(NN);

% Basis recombination

B = [zeros(1,NN-1);eye(NN-1,NN-1)];
for jj = 1:NN-1
     B(1,jj) = (-(jj+1/2)+1i*conj(lam))/(0.5-1i*conj(lam));
end


[Q,~] = qr(B,"econ");

L = L*Q;
NN = NN - 1;
L = L(:,1:NN);


if ~isempty(X)
    X=X(:);
    LL = zeros(length(X),N+1);
    LL(:,1) = exp(-X/2);
    LL(:,2) = (1 - X).*exp(-X/2);
    for jj=2:N
        k = jj - 1;
        LL(:,jj+1) = (2*k+1 - X)./(k+1).*LL(:,jj)-k/(k+1)*LL(:,jj-1);
    end
    Z = LL*Q;
else
    Z = [];
end

end







