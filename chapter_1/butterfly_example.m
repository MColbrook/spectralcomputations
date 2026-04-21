clear
close all
addpath('data')

load('butterfly_data.mat') % code for constructing the tilings and operators can be found in chapter 3: see approx_state_penrose.m
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
