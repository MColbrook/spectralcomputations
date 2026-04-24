clear
close all

% Add utils and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'utils');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'data_online');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'data_chapter7');
addpath(genpath(targetPath));

%%%%%%%%% CODE FOR PENROSE TILINGS EXAMPLE (CASE STUDY 3) %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot the generated tiles

k = 1; % level of approximation
TYPE = 1; % tile type

if TYPE == 1
    [~, loc, edges, indx_bdy] = circ_lap_dual(k-1);
    for jj=1:length(edges)
        plot(real(loc(edges(jj,:))),imag(loc(edges(jj,:))),'-','linewidth',0.5,'color',[1,1,1]*0.5);
        hold on
    end
    [L, type, loc, ang, ~] = circ_lap_tiling(k-1);
    plot_tiling(type,loc,ang,k,L)
    axis equal
    axis off
elseif TYPE == 2
    [~, loc, edges, indx_bdy] = circ_lap_dual(k-1);
    plot(real(loc),imag(loc),'.k','markersize',22-k)
    hold on
    for jj=1:length(edges)
        plot(real(loc(edges(jj,:))),imag(loc(edges(jj,:))),'-k','linewidth',1);
        hold on
    end
    axis equal
    axis off
else
    [L, loc, edges, indx_bdy] = circ_lap_dual_rhomb(k-1);
    plot(real(loc),imag(loc),'.k','markersize',22-k)
    hold on
    for jj=1:length(edges)
        plot(real(loc(edges(jj,:))),imag(loc(edges(jj,:))),'-k','linewidth',1);
        hold on
    end
    axis equal
    axis off
end


%% Approximate the spectrum
clear

k = 8; % level of approximation
TYPE = 3; % tile type
N = 5000; % truncation size

% build the matrix

if TYPE == 1
    [L, ~, loc, ~, indx_bdy] = circ_lap_tiling(k);
    X = 0:0.01:7;
elseif TYPE == 2
    [L, loc, ~, indx_bdy] = circ_lap_dual(k);
    X = 0:0.01:13;
else
    [L, loc, ~, indx_bdy] = circ_lap_dual_rhomb(k);
    X = 0:0.01:10;
end

[~,I]=sort_mi(loc);
Id = speye(length(I));
pt(I) = 1:length(I);
indx_bdy = pt(indx_bdy);
L = L(I,I);
L(:,indx_bdy)=[];
Id(:,indx_bdy)=[]; 

L = L(:,1:N);
Id = Id(:,1:N);
I = find(sum(abs(L)+abs(Id),2)<0.1);
L(I,:) = [];
Id(I,:) = [];

figure
spy(L)

E = eig(L(1:N,1:N)); % finite section

% apply CompSpec
DIST = 0*X; spec = 0*X;
pf = parfor_progress(length(X));
pfcleanup = onCleanup(@() delete(pf));
for jj = 1:length(X)
    warning('off','all')
    DIST(jj) = svds(L-X(jj)*Id,1,'smallest');
    warning('on','all')
    parfor_progress(pf);
end
for jj = 1:length(X)    
   if DIST(jj)<=0.5
       x = find (abs(X-X(jj))<DIST(jj));
       d=DIST(x);
       spec(x(d==min(d(:))))=1;
   end
end


%% Plot the results

dseig = unique(round(E'*100)/100);   
eds = ones(size(dseig));
d1   =  kron(dseig,[1 1 NaN]) + kron(eds,[1.1i 1.4i NaN]+0.1i);

dseig = X(spec==1);   
eds = ones(size(dseig));
d2   = kron(dseig,[1 1 NaN]) + kron(eds,[1.1i 1.4i NaN]+0.1i);
lw = 84*min(diff(dseig))/100;
             
figure;
tiledlayout(3,1,"TileSpacing","compact")
nexttile
plot(real(d1),imag(d1),'-','linewidth',lw, 'color','k');
set(gca,'ytick',[])
ylim([1.1-.05,1.4+.05]+0.1)
xlim([0,max(E)]);
ax = gca; ax.FontSize = 12;
title('Finite Section','interpreter','latex','fontsize',12) % Heavy spectral pollution!

nexttile
plot(X,DIST,'k')
xlim([0,max(E)]);
title('Bound on Distance to the Spectrum','interpreter','latex','fontsize',12)
ax = gca; ax.FontSize = 12;

nexttile
plot(real(d2),imag(d2),'-','linewidth',lw, 'color','k');
set(gca,'ytick',[])
ylim([1.1-.05,1.4+.05]+0.1)
xlim([0,max(E)]);
ax = gca; ax.FontSize = 12;
title('CompSpec','interpreter','latex','fontsize',12)

%% Spectral gaps - here we use precomputed high resolution approximations of the spectrum

clear
TYPE = 2;

if TYPE == 1
    load('PenTile_spec_high_res.mat');
elseif TYPE == 2
    load('PenVertex_spec_high_res.mat');
elseif TYPE == 3
    load('PenRhom_spec_high_res.mat');
end

cc = 10.^(0:-0.002:-4);
cc = cc(cc>max(DIST(spec==1)));
cc = cc(cc<max(DIST));

ng = zeros(length(cc),1);
ng2 = ng;
ng3 = ng;

ct = 1;
for cut_off=cc
    [~,I2a] = findpeaks(DIST);
    I2a = intersect(find(DIST>cut_off),I2a);
    I2a = sort(I2a);
    xx1=X(I2a);
    if min(xx1(2:end)-xx1(1:end-1))>cut_off
        ng(ct)=length(xx1);
    elseif length(xx1)==1
        ng(ct)=length(xx1);
    else
        ng(ct)=NaN;
    end

    [~,I2b] = findpeaks(FS_DIST);
    I2b = intersect(find(FS_DIST>cut_off),I2b);
    I2b = sort(I2b);
    xx2=X(I2b);
    if min(xx2(2:end)-xx2(1:end-1))>cut_off
        ng2(ct)=length(xx2);
    elseif length(xx2)==1
        ng2(ct)=length(xx2);
    else
        ng2(ct)=NaN;
    end

    [~,I2c] = findpeaks(FT_DIST);
    I2c = intersect(find(FT_DIST>cut_off),I2c);
    I2c = sort(I2c);
    xx3=X(I2c);
    if min(xx3(2:end)-xx3(1:end-1))>cut_off
        ng3(ct)=length(xx3);
    elseif length(xx3)==1
        ng3(ct)=length(xx3);
    else
        ng3(ct)=NaN;
    end
    ct = ct+1;
end

CC = gray(4); CC(end,:) = [];
CC(2:3,:) = 1.1*CC(2:3,:);
CC = flipud(CC);

figure
loglog(cc,ng,'linewidth',2,'color',CC(3,:))
hold on
ng2(isnan(ng2))=0.001;
loglog(cc,ng2,'linewidth',1,'color',CC(2,:))
loglog(cc,ng3,'linewidth',1,'color',CC(1,:))
if TYPE == 1
    title('$\#$ Gaps of $\mathrm{Sp}(H_1)$ of Length $\geq2\delta$','interpreter','latex','fontsize',18)
elseif TYPE == 2
    title('$\#$ Gaps of $\mathrm{Sp}(H_2)$ of Length $\geq2\delta$','interpreter','latex','fontsize',18)
elseif TYPE == 3
    title('$\#$ Gaps of $\mathrm{Sp}(H_3)$ of Length $\geq2\delta$','interpreter','latex','fontsize',18)
end
xlabel('$\delta$','interpreter','latex','fontsize',18)
r=max(DIST(spec==1));
plot([r,r],[1,max(max(ng2),200)],'--k','linewidth',1)
plot([1,1]*max(cc),[1,max(max(ng2),200)],'--k','linewidth',1)
legend({'$\Phi_n$ (converges)','Finite Section','Finite Tiling'},'fontsize',16,'interpreter','latex','location','southwest')
if TYPE == 1
    xlim([0.001,.1])
    ylim([1,max(max(ng2),100)])
    ax = gca; ax.FontSize = 18;
elseif TYPE == 2
    xlim([0.005,0.7])
    ylim([1,100])
    ax = gca; ax.FontSize = 18;
elseif TYPE == 3
    xlim([0.001,0.5])
    ylim([1,max(max(ng3),100)])
    ax = gca; ax.FontSize = 18;
end

%% Lebesgue measure

clear
load('PenTile_spec_high_res.mat')
r = max(max(DIST(spec==1)),0.0001);
d1 = r*(1.01.^(0:1000));
d1 = d1(d1<0.5);

LEB1 = zeros(size(d1,1),1);
for jj = 1:length(d1)
    LEB1(jj) = LEB_comp(X(1:end),DIST(1:end),d1(jj),min(X(spec==1)),max(X(spec==1)));
end
L1 = max(X(spec==1));

load('PenVertex_spec_high_res.mat')
r = max(max(DIST(spec==1)),0.0001);
d2 = r*(1.01.^(0:1000));
d2 = d2(d2<0.5);

LEB2 = zeros(size(d2,1),1);
for jj = 1:length(d2)
    LEB2(jj) = LEB_comp(X,DIST,d2(jj),min(X(spec==1)),max(X(spec==1)));
end
L2 = max(X(spec==1));


load('PenRhom_spec_high_res.mat')
r = max(max(DIST(spec==1)),0.0001);
d3 = r*(1.01.^(0:1000));
d3 = d3(d3<0.5);

LEB3 = zeros(size(d3,1),1);
for jj = 1:length(d3)
    LEB3(jj) = LEB_comp(X,DIST,d3(jj),min(X(spec==1)),max(X(spec==1)));
end
L3 = max(X(spec==1));


CC = gray(4); CC(end,:) = [];
CC = flipud(CC);

figure
semilogx(d1,LEB1/L1,'-','linewidth',4,'color',CC(1,:))
hold on
loglog(d2,LEB2/L2,'-','linewidth',1,'color',CC(2,:))
loglog(d3,LEB3/L3,'-','linewidth',2,'color',CC(3,:))

ax = gca; ax.FontSize = 18;
title('Relative Lebesgue Measure of Cover','interpreter','latex','fontsize',18)
xlabel('$\delta$','interpreter','latex','fontsize',18)
ax = gca; ax.FontSize = 18;

legend({'$T_1$','$T_2$','$T_3$'},'fontsize',16,'interpreter','latex','location','southeast')
grid minor
ax = gca; ax.FontSize = 18;


%% Box-counting dimension
clear

SamS = 10; % number of samples to avoid gridding (increase for smoother results)

for TYPE = 1:3

    if TYPE == 1
        load('PenTile_spec_high_res.mat');
    elseif TYPE == 2
        load('PenVertex_spec_high_res.mat');
    elseif TYPE == 3
        load('PenRhom_spec_high_res.mat');
    end

    r = max(DIST(spec==1))*3;
    S = cover_comp(X,DIST,r,min(X(spec==1)),max(X(spec==1)));
    d = r*(1.001.^(0:10000));
    d = d(d<0.5); d = d(:);

    dim = zeros(SamS,length(d));

    for j=1:SamS
        dim(j,:) = box_dim(S,d,0);
    end
    dim = mean(dim,1);
    
    if TYPE == 1
        d1 = d;
        dim1 = dim;
    elseif TYPE == 2
        d2 = d;
        dim2 = dim;
    else
        d3 = d;
        dim3 = dim;
    end

    dim = box_dim(S,d(d<0.05));
end


CC = gray(4); CC(end,:) = [];
CC = flipud(CC);

figure
loglog(d1,exp(log(1./d1(:)).*dim1(:)),'-','linewidth',4,'color',CC(1,:))
hold on
loglog(d2,exp(log(1./d2(:)).*dim2(:)),'-','linewidth',1,'color',CC(2,:))
loglog(d3,exp(log(1./d3(:)).*dim3(:)),'-','linewidth',2,'color',CC(3,:))

title('$N_\delta(\mathrm{Sp}(H_k))$','interpreter','latex','fontsize',18)
xlabel('$\delta$','interpreter','latex','fontsize',18)
ax = gca; ax.FontSize = 18;

legend({'$k=1$','$k=2$','$k=3$','$\sim 1/r$'},'fontsize',16,'interpreter','latex','location','southwest')
grid minor
ax = gca; ax.FontSize = 18;




























function plot_tiling(type,loc,ang,n,L)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting: shape set-up

 phi = (1+sqrt(5))/2;

 T1 = [0; 1i; exp(7i*pi/10)/phi; 0];                % B
 T2 = [0; 1i; exp(3i*pi/10)/phi; 0];                % P
 T3 = [0; exp(-1i*pi/10); exp(1i*pi/10); 0]/phi;    % W
 T4 = [0; exp(-1i*pi/10); exp(1i*pi/10); 0]/phi;    % Y

% plotting: draw shapes



    colvec = [   238 204 102; % light yellow
                 102 153 204; % light blue
                 153  68  85; % dark red
                   0  68 136; % dark blue 
             ]/255;

 tile = zeros(4,length(type));
 for j=1:length(type)
    col = colvec(type(j),:);
    switch type(j)
       case 1
         aa = T1;
       case 2
         aa = T2;
       case 3
         aa = T3;
       case 4
         aa = T4;
    end
 
    tile(:,j) = (loc(j)+ exp(1i*ang(j))*aa)*exp(1i*pi/10);
    whos tile
    % f = fill(real(tile),imag(tile),col);
    if norm(mean(tile(1:3,j)))>0.00001
        plot(mean(real(tile(1:3,j))),mean(imag(tile(1:3,j))),'k.','markersize',22-n);
    end
    % set(f,'linewidth',.01)
    hold on
 end
    tile = tile(1:3,:);
    t1 = real(mean(tile,1));
    t2 = imag(mean(tile,1));

 for j = 1:length(type)
     for k = 1:length(type)
         if L(j,k)<0
             plot(t1([j,k]),t2([j,k]),'-k','linewidth',1);
         end
     end
 end
 axis equal, axis off
end


function LEB = LEB_comp(X,dist,r,A,B)
dist = dist(:);
X = X(:);
dist(X<0)=[];
X(X<0)=[];

I = floor(r/(X(2)-X(1)));


for kk = 1: I-1

    X2 = X(kk:I:end);
    dist2 = dist(kk:I:end);
    
    S = zeros(length(X2),2);
    S(:,1) = X2 - dist2;
    S(:,2) = X2 + dist2;
    S = consolidate_int(S,10^(-14));
    % now compute the complement
    S2 = [A-10^(-12),S(1,1)+10^(-12); S(1:end-1,2)-10^(-12), S(2:end,1)+10^(-12); S(end,2)-10^(-12),B];
    S = consolidate_int(S2,10^(-14));
    LEB(kk) = sum(S(:,2)-S(:,1));
end
LEB = mean(LEB);

end


function S = cover_comp(X,dist,r,A,B)


dist = dist(:);
X = X(:);
dist(X<0)=[];
X(X<0)=[];

I = floor(r/(X(2)-X(1)));
X = X(1:I:end);
dist = dist(1:I:end);

S = zeros(length(X),2);
S(:,1) = X - dist;
S(:,2) = X + dist;

S = consolidate_int(S,10^(-14));
% now compute the complement
S2 = [A-10^(-12),S(1,1)+10^(-12); S(1:end-1,2)-10^(-12), S(2:end,1)+10^(-12); S(end,2)-10^(-12),B];
S = consolidate_int(S2,10^(-14));
end

