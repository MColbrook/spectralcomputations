clear
close all

N = 10000;
epsilon = 0.001;
k_j = -1;

%Set parameters for ODE solves   
L=10; %Scale parameter
X=1-10.^(0:(-0.001):-3); % adapt this grid around approximate eigenvalues for improved precision
X=X(:);

%Potential V(r) for Dirac op.
V=@(r) 0*r-0.8; %Modified Coulomb (actual potential is V/r)

%spectral density computation w.r.t.
f=@(r) sqrt(2)*r.*exp(-r);
g=@(r) sqrt(2)*r.*exp(-r);


%% //Discretize (*), reformulated on [-1,1], with ultraS spectral method.//

%Derivative matrix
D=ultraS.diffmat(N,1);

%Conversion matrix
S01=ultraS.convertmat(N,0,0);

%Variable coeffs and mult. matrices
x=chebfun('x'); V_x=chebfun(@(x) (1-x).*V(L*(1+x)./(1-x)));
M1=ultraS.multmat(N,L*(1+x)+V_x,0); M2=ultraS.multmat(N,(1+x).*(1-x).^2,1);
M3=ultraS.multmat(N,-L*(1+x)+V_x,0); M4=ultraS.multmat(N,1-x,0);


%Differential operators (***)
D11=S01*M1; D12=-0.5*M2*D+k_j*S01*M4; D21=0.5*M2*D+k_j*S01*M4; D22=S01*M3;
H=[D11 D12; D21 D22];

%Variable coeff for shifts
S_M=ultraS.multmat(N,L*(1+x),0); S=[S01*S_M sparse(N,N); sparse(N,N) S01*S_M];

%Right hand side
[ccn,~]=chebpts(N,1);
f_x=chebfun(@(x) f(L*(1+x)./(1-x))); g_x=chebfun(@(x) g(L*(1+x)./(1-x)));
F_vals=f_x(ccn); G_vals=g_x(ccn);
F_coeffs=chebtech1.vals2coeffs(F_vals); G_coeffs=chebtech1.vals2coeffs(G_vals);
rhs=[F_coeffs; G_coeffs];
rhs=S*rhs;


%% //Compute spectral density on grid X.//

%functions, nodes, weights for C-C quadrature approximating inner product
[ccn,ccw]=chebpts(N,1);
f_vals=2*L*f_x(ccn)./(1-ccn).^2;
g_vals=2*L*g_x(ccn)./(1-ccn).^2;

%Interlace banded blocks to get banded matrix
P=1:N; P=[P(:),P(:)+N]; P=P'; P=P(:);
PP=P*0; PP(P)=1:length(P);
H=H(P,P); S=S(P,P); rhs=rhs(P);

rho=zeros(1,length(X)); %density function
pf = parfor_progress(length(X));
pfcleanup = onCleanup(@() delete(pf));

for n=1:length(X)
	%Step 1: apply resolvent at node z=X(n)+eps*1i
	z=X(n)+1i*eps;
	u_coeffs=(H-z*S)\rhs;
    u_coeffs=u_coeffs(PP);

	%solution values from coeffs
	u_vals=chebtech1.coeffs2vals(u_coeffs(1:N));
	w_vals=chebtech1.coeffs2vals(u_coeffs((N+1):end));

	%Step 2: evaluate inner product with RHS via  Clenshaw-Curtis rule
	dot=ccw*(conj(f_vals).*u_vals)+ccw*(conj(g_vals).*w_vals);
	rho(n)=imag(dot)*epsilon;
	parfor_progress(pf);
end


%% Plot measure
figure
loglog(1-X,rho)

xlabel('$1-x$','interpreter','latex','fontsize',18)
title('$\upsilon_v^\epsilon(x)$','interpreter','latex','fontsize',18)

ax=gca; ax.FontSize=17;
