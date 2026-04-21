clear
close all

addpath(genpath('./routines'))

%% Construct the matrix

phi = 0.69; % not 0.59 due to scaling used in Hamiltonian construction
[H,f,Pos] = ham_construct(phi,117998); % slowest part of the code!

figure % plot the function f
Pos=Pos(:); % position of lattice points in real space
plot(max(f-transpose(1:length(f)),0))
pause(0.1)

%% Plot the approximate eigenvectors

Z = 3.04294150;
Z = 1.395230; % spectral parameter where we hunt for approximate state

n = 100000;
H1=(H(1:f(n),1:n)-Z*speye(f(n),n));
[~,~,psi1] = svds(H1,1,'smallest');
En = H(1:f(n),1:n)*psi1;
Z = En(1:n)'*psi1; % Rayleigh quotient

figure
scatter(real(Pos(1:length(psi1))),imag(Pos(1:length(psi1))),9,0*(abs(psi1)),'MarkerEdgeColor',[1,1,1]*0,'Linewidth',0.1);
hold on
scatter(real(Pos(1:length(psi1))),imag(Pos(1:length(psi1))),8,(abs(psi1)),'filled');
axis equal; axis([-20,20,-20,20]); ax = gca; ax.YTick=[]; ax.XTick=[]; set(gca, 'Layer', 'Top'); box on
colormap(flipud(bone))
clim([0*max((abs(psi1))),0.7*max((abs(psi1)))])

norm(H(1:f(n),1:n)*psi1-Z*speye(f(n),n)*psi1) % error bound
