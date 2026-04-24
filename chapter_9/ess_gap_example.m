clear
close all

%%%%%%%%%%%%%%%%%%%%%%% CODE FOR TABLE 9.2 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up the matrix discretisation of the operator A0

M = 50; % size of truncation is 2*M+1
N = 2*(M+10)+1; % padding to allow easy construction using powers of shift
S = speye(N);
S = [sparse(N,1),S;sparse(1,N+1)];
S = S(1:N,1:N); % shift operator

A = (S*S*S*S+S')/2;
ORD = -(M+10):(M+10);
A = A+  spdiags(repmat((3i./(abs(ORD)+1))',1,1),1,N,N); % add perturbation

%% Compute eigenvalues

[V,E] = eig(full(A(11:end-10,11:end-10)),'vector');
R = sort(abs(E),'descend');

R(1:5)


