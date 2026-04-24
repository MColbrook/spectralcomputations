clear
close all


%% Build the discretised operator
n1 = 300;
n2 = 150;

dx = spdiags(sqrt((0:500)'/2),1,500,500);
mx = dx+dx';
dx = dx-dx';
A = -dx*dx +1i*mx;

A = A((n2+1):(n1+n2),(n2+1):(n1+n2));

%% Spectral data
p1 = Num_Range(A,200); % numerical range using Johnson's method
E = sqrt(1i)*(2*(0:100)+1); % analytic eigenvalues


%% Plot the results

T1 = sprintf('We_airy_%d_%d.pdf',n1,n2);
T2 = sprintf('$n_1=%d$, $n_2=%d$',n1,n2);


figure
fill(real([p1(:);p1(1)]),imag([p1(:);p1(1)]),[1,1,1]*0.65)
axis equal
axis([-10,50,-20,20])

ax=gca; ax.FontSize=18;
title(T2,'interpreter','latex','fontsize',26)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',26)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',26)

exportgraphics(gcf,T1,'ContentType','vector','BackgroundColor','none','Resolution',300,'Colorspace','gray')

function [p] = Num_Range(A,k)
theta=(0:k-1)*2*pi/k;
p=0*theta;

for j=1:length(theta)
    H=exp(1i*theta(j))*A;
    H=(H+H')/2;
    [V,E]=eig(full(H));
    I=find(diag(E)==max(diag(E)));
    p(j)=(V(:,I(1)))'*(A*V(:,I(1)));
end

end