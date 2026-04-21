clear
close all
addpath(genpath('FEM_data'));
load('3D_maxwell.mat') % FEM matrices computed using DOLFIN

v=vecnorm(B_edge);
I = find(v>10^(-12));
B_edge = B_edge(I,I); B_edge = (B_edge+B_edge')/2;
A_edge = A_edge(I,I); A_edge = (A_edge+A_edge')/2;
A0_edge = A0_edge(I,I)*1i; A0_edge = (A0_edge+A0_edge')/2;

%% Zimmermann and Mertins trick (may need extra memory depending on machne)

tup = 0.05;
tlow = 2.8;

x = tup;
M1 = A_edge-x*(A0_edge'+A0_edge)+abs(x)^2*B_edge;
M2 = A0_edge-x*B_edge;

S1 = eigs(M2,M1,16,2,'Tolerance',1e-7);

x = tlow;
M1 = A_edge-x*(A0_edge'+A0_edge)+abs(x)^2*B_edge;
M2 = A0_edge-x*B_edge;

S2 = eigs(M2,M1,16,-20,'Tolerance',1e-12);

b_bound = sort(tup+1./S1);
a_bound = sort(tlow+1./S2);

%% Without resolvent

E = [1.14,1.54,2.08,2.23,2.32,2.34,2.4,2.6];
EE =E;
mult = [1,2,3,2,1,2,1,3];

pf = parfor_progress(length(EE));
pfcleanup = onCleanup(@() delete(pf));

for jj = 1:length(EE)
    x = EE(jj);
    M1 = A_edge-x*(A0_edge'+A0_edge)+abs(x)^2*B_edge;
    M2 = B_edge;
    [V,~] = eigs(M1,M2,mult(jj),0,'Tolerance',1e-12);
    RHO = real(dot(V,A0_edge*V));
    EP = real(dot(V,A_edge*V-(A0_edge'+A0_edge)*V.*RHO+B_edge*V.*RHO.^2));
    Ecomp_edge{jj} = RHO;
    Err_edge{jj} = EP;
	parfor_progress(pf);
end

%% Table of data

format long

E1=[Ecomp_edge{:}]';
E2=sqrt([Err_edge{:}]');
[~,I]=sort(E1);

E1=E1(I);
E2=E2(I);
c=[round(10^6*E1)/10^6,ceil(E2*1000000)/1000000,floor(10^6*a_bound((a_bound>1)&(a_bound<2.7)))/10^6,ceil(10^6*b_bound((b_bound>1)&(b_bound<2.7)))/10^6,...
    ceil(10^6*(b_bound((b_bound>1)&(b_bound<2.7))-a_bound((a_bound>1)&(a_bound<2.7)))/2)/10^6];

fprintf('%9f & %9f & %9f & %9f & %9f \\\\ \n', c')





