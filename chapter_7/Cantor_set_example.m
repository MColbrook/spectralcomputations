clear
close all

addpath(genpath('./routines'))

%%%%%%%%%%%%%%%%% CODE FOR THE CANTOR SET EXAMPLES %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot the X_n that generate the Cantor set

N = 5; % number of levels in the plot
X = cell(N+2,1); % cell array of the intervals

X{1} = [0,1]; % n = 0
for n = 1:N
    X{n+1} = [X{n};X{n}+2]/3; % generate the sets X_n
end

f = figure;
plot(real(X{1}),0*X{1},'-','linewidth',20, 'color','k');
hold on
for n = 1:N
    plot(real(X{n+1}'),0*X{n+1}'-n/2,'-','linewidth',20, 'color','k'); % plot the unions of intervals X_n
end
set(gca,'ytick',fliplr(-(0:N)/2),'yticklabels',["$n=6$","$n=5$","$n=4$","$n=3$","$n=2$","$n=1$","$n=0$"])
set(gca,'YColor',[1 1 1])
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';
yaxisproperties.TickLabelColor = 'k';
title('$X_n$','interpreter','latex','fontsize',18)
ylim([-(N)/2-0.25,0.2])
xlim([-0.002,1])
ax=gca; ax.FontSize=18;
f.Position = [360.0000   97.6667  812.3333  420.0000];
box off

%% Approximate the box-counting dimension - naive without proper scaling fails

N = 25; % number of levels in the plot
X = cell(N+2,1); % cell array of the intervals

X{1} = [0,1]; % n = 0

for n = 1:N
    X{n+1} = [X{n};X{n}+2]/3; % generate the sets X_n
end

delta = 2.^(0:-0.2:-50);

dim = box_dim(X{5+1},delta,0);
NC1 = exp(log(1./delta(:)).*dim(:));

dim = box_dim(X{10+1},delta,0);
NC2 = exp(log(1./delta(:)).*dim(:));

dim = box_dim(X{15+1},delta,0);
NC3 = exp(log(1./delta(:)).*dim(:));

% plot the results
CC = gray(4); CC(end,:) = [];
CC = flipud(CC);

figure
loglog(delta,NC1,'-','linewidth',4,'color',CC(1,:))
hold on
loglog(delta,NC2,'-','linewidth',1,'color',CC(2,:))
loglog(delta,NC3,'-','linewidth',2,'color',CC(3,:))
loglog(delta,(1./delta).^(log(2)/log(3)),'k--','linewidth',1)
loglog(delta(delta<10^(-5)),1./delta(delta<10^(-5)),'k:','linewidth',2)
text(10^(-9),10.^9*4,'wrong scaling','interpreter','latex','fontsize',16,'Rotation',-40)
text(10^(-5),10.^3/2,'correct scaling','interpreter','latex','fontsize',16,'Rotation',-28)
xlim([10^(-12),1])

title('$N_\delta(X_k)$','interpreter','latex','fontsize',18)
xlabel('$\delta$','interpreter','latex','fontsize',18)
ax = gca; ax.FontSize = 18;

legend({'$k=5$','$k=10$','$k=15$','$(\frac{1}{\delta})^{\frac{\log(2)}{\log(3)}}$','${1}/{\delta}$'},'fontsize',16,'interpreter','latex','location','southwest','NumColumns',1)
grid minor
ylim([1,10^10])
ax = gca; ax.FontSize = 18;

%% Approximate the box-counting dimension - adaptive choice (in this simple example, just use large enough N)

dim = box_dim(X{N+1},delta,0);
NC = exp(log(1./delta(:)).*dim(:));

figure
loglog(delta,NC,'k-','linewidth',3)
hold on
loglog(delta,(1./delta).^(log(2)/log(3)),'k--','linewidth',1)
xlim([10^(-12),1])

title('$N_\delta(X)$','interpreter','latex','fontsize',18)
xlabel('$\delta$','interpreter','latex','fontsize',18)
ax = gca; ax.FontSize = 18;

grid minor
ylim([1,10^10])
ax = gca; ax.FontSize = 18;

%% Hausdorff dimension plot

N = 19; % number of levels in the plot - increase for more of the plot
X = cell(N+1,1); % cell array of the intervals

X{1} = [0,1]; % n = 0

for n = 1:N
    X{n+1} = [X{n};X{n}+2]/3; % generate the sets X_n
end

er0 = (1/3)^N;
n1max = min(ceil(log2(1/er0)),30);
dimH = Haus_dim(X{N+1},5:n1max,5:5:20,0,1.3);

CC = gray(5); CC(end,:) = [];
CC = flipud(CC);
figure
p1 = plot(5:30,dimH(:,1),'-^','markersize',7,'color',CC(1,:),'linewidth',1);
p1.MarkerEdgeColor = [0 0 0];
p1.MarkerFaceColor = CC(1,:);
hold on
p2 = plot(5:30,dimH(:,2),'-s','markersize',7,'color',CC(2,:),'linewidth',1);
p3 = plot(5:30,dimH(:,3),'-d','markersize',7,'color',CC(3,:),'linewidth',1);
p3 = plot(5:30,dimH(:,4),'-o','markersize',7,'color',CC(3,:),'linewidth',1);
p2.MarkerEdgeColor = [0 0 0];
p2.MarkerFaceColor = CC(2,:);
p3.MarkerEdgeColor = [0 0 0];
p3.MarkerFaceColor = CC(3,:);
ax=gca; ax.FontSize=18;
xlabel('$n_1$','interpreter','latex')
plot([0,30],[0,0]+log(2)/log(3),':k','linewidth',1)
ylim([0.55,0.75])
xlim([5,30])
title('$\Gamma_{n_2,n_1}$','interpreter','latex')

hleg = legend({'$n_2=5$','$n_2=10$','$n_2=15$','$n_2=20$'},...
    'interpreter','latex', 'fontsize',14,'location','best');
grid minor


