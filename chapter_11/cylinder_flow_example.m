clear
close all

% Add Koopman algorithms and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'Koopman_algorithms');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'Koopman_datasets');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%%% CODE FOR CYLINDER FLOW EXAMPLE %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Vorticity_data.mat')

%% Perform EDMD with POD modes (adjoint of exact DMD!)
r = 47;
M = 24*5;
ind = (1:M);
X = VORT(:,ind);
Y = VORT(:,ind+1);
[U,S,~] = svd(X,'econ');
r = min(rank(S),r);
U = U(:,1:r);

PX = X'*U;
PY = Y'*U;
K = PX\PY;
[V,LAM] = eig(K,'vector');

Phi = transpose((PX*V)\(X')); % Koopman modes
[~,I] = sort(abs(1-LAM),'ascend'); % reorder modes
Phi = Phi(:,I); LAM = LAM(I);

%% Plot the eigenvalues
figure
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-k')
hold on
plot(real(LAM),imag(LAM),'k.','markersize',18)
axis equal
axis([-1.15,1.15,-1.15,1.15])

for j=0:5
    text(1.05*real(LAM(max(1,2*j))),1.05*imag(LAM(max(1,2*j))),sprintf('%d',j),'interpreter','latex','fontsize',13)

end
text(1.05*real(LAM(max(1,2*7)))+0.03,1.05*imag(LAM(max(1,2*7)))-0.03,'$\ddots$','interpreter','latex','fontsize',13,'rotation',-5)

title('EDMD Eigenvalues','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

%% Predictive power
b = Phi\VORT(:,1);
Er = zeros(1000,1);
m = mean(VORT(:,1:M),2);

for j=1:1000
    Er(j)=norm(Phi*(LAM.^(j).*b)-VORT(:,j+1))/norm(VORT(:,j+1)-m);
end

figure
loglog(Er,'k','linewidth',2)
hold on
plot([120,120],[10^(-9),max(Er)],'--k','linewidth',2)
title('Relative Prediction Error','interpreter','latex','fontsize',18)
xlabel('Number of time steps ($n$)','interpreter','latex','fontsize',18)
xx = [0.5 0.6]+0.05;    yy = [0.8 0.7];
annotation('textarrow',xx,yy,'String','Extent of snapshot data','interpreter','latex','fontsize',16,'Linewidth',1)
ax=gca; ax.FontSize=18;

%% Plot modes
ct = 0;
for j=[1,3,5,7]
    figure
    
    C = zeros(800*200,1)+NaN;
    C(II) = real(Phi(:,j));
    
    C = reshape(C,[800,200]);
    
    if mod(floor((j)/2),2)==1
        C=(C+fliplr(C))/2;
    else
        C=(C-fliplr(C))/2;
    end
    vv=-0.025:0.05:1;
    
    a = prctile(C(~isnan(C)),100);
    b = prctile(C(~isnan(C)),0);
    a = max(abs(a),abs(b)); b = -a;
    C(~isnan(C)) = max(C(~isnan(C)),b);
    C(~isnan(C)) = min(C(~isnan(C)),a);

    c=(a+b)/2;
    
    [~,h]=contourf(Xgrid,Ygrid+0.06,C,vv*(a-b)+b,'edgecolor','k');
    h.LineWidth = 0.2;
    colormap(brighten(gray,0.5))
    colorbar
    clim([b,a])

    axis equal
    hold on
    fill(1.1/2*cos(0:0.01:2*pi),1.1/2*sin(0:0.01:2*pi),'w','edgecolor','none')
    plot(1.1/2*cos(0:0.01:2*pi),1.1/2*sin(0:0.01:2*pi),'k','linewidth',1)
    xlim([-2,10])
    title(sprintf('Mode %d',ct),'interpreter','latex','fontsize',16)
    ylim([-2,2])
    xlabel('$x/D$','interpreter','latex','fontsize',18)
    ylabel('$y/D$','interpreter','latex','fontsize',18)
    hold off
    ax=gca; ax.FontSize=18;
    ct = ct +1;
end


