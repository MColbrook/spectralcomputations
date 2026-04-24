clear
close all

%%%%%%% CODE FOR THE PERTURBED LAURENT EXAMPLE (FIGURE 8.7) %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up the matrix discretisation of the operator A0

M = 100;
N = 2*(M+10)+1; % padding to allow easy construction using powers of shift
S = speye(N);
S = [sparse(N,1),S;sparse(1,N+1)];
S = S(1:N,1:N); % shift operator

A = (S*S*S*S+S')/2;

p2 = Num_Range(A(11:end-10,11:end-10),200); % numerical range of A0 using Johnson's method

%% Add the perturbation and compute eigenvalues

ORD = -(M+10):(M+10);
A = A+spdiags(repmat((3i./(abs(ORD)+1))',1,1),1,N,N);

[V,E] = eig(full(A(11:end-10,11:end-10)),'vector');

II = speye(size(A));
res = vecnorm(A(:,11:end-10)*V-II(:,11:end-10)*V*diag(E));

A = A(11:end-10,11:end-10);
p1 = Num_Range(A,200); % numerical range of A using Johnson's method

%% Plot the results

figure
fill(real([p1(:);p1(1)]),imag([p1(:);p1(1)]),[1,1,1]*0.9,'LineStyle','none')
hold on
fill(real([p2(:);p2(1)]),imag([p2(:);p2(1)]),[1,1,1]*0.65,'LineStyle','none')

theta = -pi:0.01:pi;
t = exp(1i*theta);
sy = (t.*t.*t.*t+conj(t))/2;
plot(real(sy),imag(sy),'k','linewidth',1) % essential spectrum (symbol curve)

plot(real(E(res>0.1)),imag(E(res>0.1)),'ok') % "spurious"
plot(real(E(res<0.1)),imag(E(res<0.1)),'xk','markersize',12,'linewidth',2) % "non-spurious"
axis equal

ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
legend({'$W(A)$','$W_e(A)$','$\mathrm{Sp}_{\mathrm{ess}}(A)$','Pollution','Reliable'},'fontsize',14,'interpreter','latex','location','northeast')
ax=gca; ax.FontSize=18;

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

