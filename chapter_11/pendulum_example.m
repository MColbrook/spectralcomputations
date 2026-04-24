clear
close all

% Add Koopman algorithms and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'Koopman_algorithms');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'Koopman_datasets');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%%%%%%%% CODE FOR FIGURE 11.29 %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up the experiment

N = 750; % number of delay embeddings + 1
M1 = 250; % number of x data points (500 for book)
M2 = 250; % number of y data points (500 for book)

g = @(x,y) exp(5i*x)./cosh(y); % observable

%% Collect the data
[PX,PY,W,x1,x2] = pendulum_matrix(M1,M2,N,5,g); % computes matrices using time-delay embedding

%% Perform riggedDMD
TH2 = -pi:0.005:pi;                     % angles for spectral measure
epsilon = 0.05;                         % smoothing parameter
order = 6;                              % order of kernel
g_coeffs = zeros(size(PX,2),1); g_coeffs(1)=1;   % cofficients of g in dictionary expansion

[~,xi] = riggedDMD(PX,PY,W,epsilon,[],[],'order',order,'g_coeffs',g_coeffs,'TH2',TH2);
xi = (max(xi,0.000001)/sum(xi*(TH2(2)-TH2(1)))); % spectral measure normalised to be a probability measure

[~,locs] = findpeaks(xi);
locs = locs(TH2(locs)>-0.005);
TH = TH2(locs); % angles for generalised eigenfunction
[gTH,~] = riggedDMD(PX,PY,W,epsilon,TH,[],'order',order,'g_coeffs',g_coeffs);

%% Plot the results
close all
cc = [1,2,5,6,9,10,12,13,14,15,16,18];
locs2 = locs(cc);

f = figure % spectral measure plot
plot(TH2,xi,'k','linewidth',2)
xlim([-0.0015,pi])
ax=gca; ax.FontSize=18;
hold on
plot(TH2(locs2),xi(locs2),'ok','markersize',12,'linewidth',1)
grid on
box on
ax=gca; ax.FontSize=18;
f.Position = [40   97.6667  560.*2  420.];
xlabel('$\theta$','interpreter','latex','fontsize',18)


for jj = 1:length(cc)
    figure % generalised eigenfunction plot
    C = PX(:,1:400)*gTH(1:400,cc(jj));
    C = C-mean(C(:));

    toPlot = real((reshape(C,length(x2),length(x1))));
    contourf(x1, x2 , toPlot ,60,'edgecolor','none');
    colormap(brighten(brewermap([],'RdYlBu'),0))
    axis equal on;   view(0,90);    ylim([-4,4]);   xlim([-pi,pi])
    clim([-max(abs(toPlot(:)))*0.8,max(abs(toPlot(:)))*0.8])
    ax=gca; ax.FontSize=16;
    set(gca,'YDir','normal')
    colormap(brighten(gray,0.5))
    pause(0.1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PSI_X,PSI_Y,W,x1,x2] = pendulum_matrix(M1,M2,N,L,g)
% computes data matrices for the pendulum
delta_t=0.4;
x1=linspace(-pi,pi,M1+1);   x1=x1+(x1(2)-x1(1))/2; x1(end)=[];
x2=linspace(-L,L,M2+1); x2=x2+(x2(2)-x2(1))/2; x2(end)=[];
[X1,X2] = meshgrid(x1,x2);  X1=X1(:); X2=X2(:);
M=length(X1); % number of data points

Y0=[X1(:)';X2(:)'];
[~,Y]=ode45(@(t,y) pensystem(y,length(X1)),[0.000001 (1:N)*delta_t],Y0);

XX = mod([X1(:)';Y(2:end,1:2:end)]+pi,2*pi)-pi;
YY = [X2(:)';Y(2:end,2:2:end)];

Y = transpose(g(XX,YY));
PSI_X = Y(:,1:end-1);
PSI_Y = Y(:,2:end);

W=zeros(M,1)+(x1(2)-x1(1))*(x2(2)-x2(1));
end



function dydt = pensystem(y,n)
y = reshape(y,[],n);
dydt = [y(2,:);-sin(y(1,:))];

% Linearize output.
dydt = dydt(:);
end