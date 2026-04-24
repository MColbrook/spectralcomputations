clear
close all

% Add Koopman algorithms and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'Koopman_algorithms');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'Koopman_datasets');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%%% CODE FOR ARCTIC SEA ICE EXAMPLE %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load the data

load('ICE_DATA.mat')
delays = 6; % number of time delays
X = DATA(:,delays:end-1);
for jj = delays-1:(-1):1
    X = [X;DATA(:,jj:end-(delays-jj+1))];
end


%% Run the algorithm

[G,K,L,PX] = kResDMD(X(:,1:end-1),X(:,2:end),'type',"Gaussian");
[V,LAM,W] = eig(K,'vector'); W = conj(W);
R = real(sqrt(dot(V,L*V+V*diag(abs(LAM)).^2-K'*V*diag(LAM)-K*V*diag(conj(LAM)))./dot(V,V))); % error bounds
[~,I] = sort(R,'ascend');
V = V(:,I); LAM = LAM(I); W = W(:,I); R = R(I);

PXr = PX*W;

%% Error bounds of EDMD eigenvalues
[~,I] = sort(R,'descend');
V = V(:,I); LAM = LAM(I); W = W(:,I); R = R(I);

figure
n2 = length(LAM);
plot(angle(LAM(end-n2+1:end)),log(abs(LAM(end-n2+1:end))),'.','markersize',18,'color',[1,1,1]*0.7);
n2 = 17;
hold on
plot(angle(LAM(end-n2+1:end)),log(abs(LAM(end-n2+1:end))),'ko','markersize',10);
plot(angle(LAM(end-n2+1:end)),log(abs(LAM(end-n2+1:end))),'.k','markersize',18)
xlabel('$\mathrm{arg}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\log(|\lambda|)$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
xlim([-pi,pi])
ylim([-0.01,0.002])
grid off
box on

%% Koopman modes

[~,I] = sort(R,'ascend');
V = V(:,I); LAM = LAM(I); W = W(:,I); R = R(I);
Phi = transpose((PX*W)\(X(1:97877,1:end-1)'));


for jj = [3,1,4,13,14,16]
    figure
    u = abs(Phi(1:size(DATA,1),jj));
    u = real(u*exp(1i*mean(angle(u))));
    
    v = zeros(432*432,1)+NaN;
    v(nLAND) = u(:);
    v = reshape(v,[432,432]);
    imagesc(v,'AlphaData',~isnan(v))
    colormap(flipud(brighten(gray,0.5)))
    set(gca,'Color',[1,1,1]*0.6)
    colorbar('southoutside')
    if round(angle(LAM(jj))/pi*6) ==1
        title(sprintf('$\\lambda=%.3fe^{\\pi\\mathrm{i}/6}$, $\\mathrm{res}^*=%.3f$',abs(LAM(jj)/max(abs(LAM))),R(jj)),'interpreter','latex','fontsize',18)
    elseif round(angle(LAM(jj))/pi*6) ==2
        title(sprintf('$\\lambda=%.3fe^{2\\pi\\mathrm{i}/6}$, $\\mathrm{res}^*=%.3f$',abs(LAM(jj)/max(abs(LAM))),R(jj)),'interpreter','latex','fontsize',18)
    else
        title(sprintf('$\\lambda=%.3f$, $\\mathrm{res}^*=%.3f$',abs(LAM(jj)/max(abs(LAM))),R(jj)),'interpreter','latex','fontsize',18)
    end
    axis equal
    axis tight
    grid off
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
    set(gca, 'XTick', [], 'YTick', [])
    pause(0.1)
end



