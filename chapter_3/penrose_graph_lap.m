clear
close all
addpath(genpath('./routines'))

k = 6; % larger k leads to larger truncation size

%% Finite section approximations

[H_fs, loc, ~, ~] = circ_lap_dual_rhomb_fs(k);
[~,I]=sort_mi(loc);
H_fs = H_fs(I,I);

d = symrcm(H_fs);
H_fs = H_fs(d,d);

E_fs = eig(H_fs);


%% Finite tiling approximations

[H_ft, loc, ~, ~] = circ_lap_dual_rhomb(k);
[~,I]=sort_mi(loc);
H_ft = H_ft(I,I);

d = symrcm(H_fs);
H_ft = H_ft(d,d);

E_ft = eig(H_ft);


%% Periodic approximations

H_per = per_penrose_construct(k);
E_per = eig(H_per);

%% CompSpec

[L, loc, ~, indx_bdy] = circ_lap_dual_rhomb(k);
[~,I]=sort_mi(loc);
Id = speye(length(I));
pt(I) = 1:length(I);
indx_bdy = pt(indx_bdy);
L = L(I,I);


figure
L(:,indx_bdy)=[]; % rectangular truncation
Id(:,indx_bdy)=[];
spy(L)

delta = 1/100;
X = -10*delta:delta:(9.2+2*delta);
DIST = 0*X;
spec = 0*X;


C2 = L'*L;
C1 = Id'*L;
C0 = Id'*Id;

p = symrcm(C2-2*15*C1+15^2*C0);
C2 = C2(p,p);
C1 = C1(p,p);
C0 = C0(p,p);


pf = parfor_progress(length(X));
pfcleanup = onCleanup(@() delete(pf));

warning ('off','all');
for jj = 1:length(X)
    DIST(jj) = sqrt(eigs(C2-2*X(jj)*C1+X(jj)^2*C0,1,'smallestabs','maxit',2000,'tol',10^(-3)));
    % DIST(jj) = svds(L-X(jj)*Id,1,'smallest','tolerance',delta/10);
    parfor_progress(pf);
end
warning ('on','all');

for jj = 1:length(X)    
    if DIST(jj)<=0.5
        x = find (abs(X-X(jj))<DIST(jj));
        d=DIST(x);
        spec(x(d==min(d(:))))=1;
    end
end
figure
plot(X,DIST)


%% Plot the results

PER_spec = uniquetol(E_per,1e-6);
pereig = PER_spec';   
eper = ones(size(pereig));

FT_spec = uniquetol(E_ft,1e-6);
FS_spec = uniquetol(E_fs,1e-6);

dseig = X(spec==1);   
eds = ones(size(dseig));

fteig = FT_spec';
eft = ones(size(fteig));

fseig = FS_spec';
efs = ones(size(fseig));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING PARAMETERS - FULL SPECTRUM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basis = 84;                    % adust to plot over [-.15, 6.15]
ax  = [0.05 0.05 .9 .5];       % to set up plotting axes
gray = .65*[1 1 1];

xw = min(diff(dseig));         % xw = 1e-4 
lw = basis*xw;                 % "true" width for lines of width xw.
xspan = 6;                     % linewidth set for width 6

zoom = [-.15 9.2; 4.4 5; -0.006 .102; -0.0008 .0104];

fsz   = kron(fseig,[1 1 NaN]) + kron(efs,[1.1i 1.4i NaN]+0.1i);
ftz   = kron(fteig,[1 1 NaN]) + kron(eft,[0.8i 1.1i NaN]+0.05i);
perz   = kron(pereig,[1 1 NaN]) + kron(eper,[0.5i 0.8i NaN]);

dsz   = kron(dseig,[1 1 NaN]) + kron(eds,[0.1i 0.4i NaN]);


for k=1:2
    f = figure;
    
    tiledlayout(4,1,"TileSpacing","compact")
    nexttile
    plot(real(fsz),imag(fsz),'-','linewidth',lw, 'color','k');
    set(gca,'ytick',[])
    ylim([1.1-.05,1.4+.05]+0.1)
    xlim(zoom(k,:));
    title('Finite Section','interpreter','latex')
    
    nexttile
    plot(real(ftz),imag(ftz),'-','linewidth',lw, 'color','k');
    set(gca,'ytick',[])
    ylim([1.1-.05,1.4+.05]+0.05-0.3)
    xlim(zoom(k,:));
    title('Finite Tiling','interpreter','latex')
    
    nexttile
    plot(real(perz),imag(perz),'-','linewidth',lw, 'color','k');
    set(gca,'ytick',[])
    ylim([.5-.05,.8+.05])
    xlim(zoom(k,:));
    title('Periodic Approximation','interpreter','latex')
    
    nexttile
    plot(real(dsz),imag(dsz),'-','linewidth',lw, 'color','k');
    set(gca,'ytick',[])
    ylim([.1-.05,.4+.05])
    xlim(zoom(k,:));
    title('\texttt{CompSpec}','interpreter','latex')
    f.Position=[360.0000   50.3333  560.0000  500];

    f = figure;
    tiledlayout(4,1,"TileSpacing","compact")
    nexttile
    plot(X,DIST, 'color','k')
    xlim(zoom(k,:));
    ylim([-0.05,.45])
    title('$\Phi_{n}(z,H_0)$','interpreter','latex')
    f.Position=[360.0000   50.3333  560.0000  590.6667];
end
