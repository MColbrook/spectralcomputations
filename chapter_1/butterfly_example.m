clear
close all
addpath('data')

%% Example of how to produce data (increase n and resolution of phi as wanted, to reduce grid size use points close to finite section eigenvalues)
n = 200;
X = -4.3:0.02:4.3;
X2 = 0:0.02:4.3;

C1 = [];
C2 = [];
for phi = 0:0.02:1
    [H,f,Pos] = ham_construct(phi*0.69/0.59,1000);
    [V,E] = eig(full(H(1:n,1:n))); E = diag(E);
    C1 = [C1;E(:) + phi*1i];
    
    DIST = 0*X2;
    for jj = 1:length(X)
        DIST(jj) = svds(H(1:f(n),1:n)-X(jj)*speye(f(n),n),1,'smallest');
    end
    DIST = [fliplr(DIST(2:end)),DIST];
    
    for jj = 1:length(X)    
        if DIST(jj)<=0.5
            x = find (abs(X-X(jj))<DIST(jj));
            d = DIST(x);
            C2 = [C2,X(x(d==min(d(:)))) + phi*1i];
        end
    end
    plot(E,0*E+phi,'.k','markersize',2)
    hold on
    pause(0.01)
end

% plot the results

f=figure
plot(real(C1),imag(C1),'.','markersize',4,'color',[1,1,1]*0.7)
ylim([0,1])
xlim([-4.3,4.3])
ax = gca; ax.FontSize = 14;
f.Position=[360.0000  345.6667  560.0000  272.3333];

f=figure
plot(real(C2),imag(C2),'.','markersize',4,'color',[1,1,1]*0.7)
ylim([0,1])
xlim([-4.3,4.3])
ax = gca; ax.FontSize = 14;
f.Position=[360.0000  345.6667  560.0000  272.3333];

%% Resolved plot
clear
load('butterfly_data.mat') % code for constructing the tilings and operators can be found in chapter 3
C3 = setdiff(round(imag(C1*400))*1i/400+round(real(C1*100))/100,round(imag(C2*400))*1i/400+round(real(C2*100))/100);
C1 = round(C2*400)/400;
C2 = round(C2*400)/400;
C2 = setdiff(C2,C3);


f=figure
plot(real(C3),imag(C3),'.k','markersize',2)
hold on
plot(real(C1),imag(C1),'.','markersize',2,'color',[1,1,1]*0.7)
ylim([0,1])
xlim([-4.3,4.3])
ax = gca; ax.FontSize = 14;
f.Position=[360.0000  345.6667  560.0000  272.3333];
exportgraphics(gcf,'CHAP1_butterfly1.pdf','BackgroundColor','none','Resolution',500,'Colorspace','gray')

f2=figure

plot(real(C2),imag(C2),'.k','markersize',1)
ylim([0,1])
xlim([-4.3,4.3])
ax = gca; ax.FontSize = 14;
f2.Position=[360.0000  345.6667  560.0000  272.3333];
exportgraphics(gcf,'CHAP1_butterfly2.pdf','BackgroundColor','none','Resolution',500,'Colorspace','gray')
