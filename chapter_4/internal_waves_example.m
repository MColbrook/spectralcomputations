clear
close all

% Add utils and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'utils');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'data_online');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%%% CODE FOR EXAMPLE IN SECTION 4.5.3 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
epsilon = 0.01;
order = 6;
X2 = 0.790.^2; % square root of this is lambda in the example

%% Set up the functions and domain for the example
hh=1;
n=70;

dom = ultraSEM.Domain.quad([0 0 ; 1+0.5 0 ; 1 1; 0 1]);
if hh>0
    dom=refine(dom,hh);
end

rng(1)
g = randnfun2(0.3,[0,1+0.5,0,1]);
rhs = @(x,y) g(x,y);  % function g
rhs2 = ultraSEM.Sol(rhs, n, dom);

%% Compute f
pdo={{@(x,y) -1+0*x+0*y,@(x,y) 0*x+0*y,@(x,y) -1+0*x+0*y},...
               {@(x,y) 0*x+0*y,@(x,y) 0*x+0*y},@(x,y) 0*x+0*y};
op=ultraSEM(dom, pdo, rhs, n); bc = 0;
f=op \ bc;


c = sqrt(sum2(f.*rhs2));
rhs = @(x,y) g(x,y)/c;  % normalise
rhs2 = (1/c)*rhs2;


%% Generalised eigenfunction

I = 0;
[pts,al]=rational_kernel(order,'equi');
for ii=1:order
    pdo={{@(x,y) X2-pts(ii)*epsilon,@(x,y) 0,@(x,y) X2-pts(ii)*epsilon-1+0},...
        {@(x,y) 0,@(x,y) 0},@(x,y) 0};
    op=ultraSEM(dom, pdo, rhs, n);
    sol=op \ 0; 
    I=I-(1/pi)*(al(ii)*sol);
end

%% Plot the generalised eigenfunction

U1 = -imag(I);
U2 = abs(diff(I,1));

figure
pcolor(U1)
view(2)
hold on
plot([0,1+0.5,1+0.5,1,1,0,0,0,1],[0,0,0,1,1,1,1,0,0],'k','linewidth',1)
clim([-norm(U1),norm(U1)])
colormap gray
axis equal tight
axis off

figure
pcolor(U2)
view(2)
hold on
plot([0,1+0.5,1+0.5,1,1,0,0,0,1],[0,0,0,1,1,1,1,0,0],'k','linewidth',1)
clim([-norm(U2),norm(U2)])
colormap gray
axis equal tight
axis off

%% Compute spectral measures

X = 0:0.01:1;
ep_vec = 0.1;
lambda = 0.1; % lengthscale of random function
dd = 0.5;

g = randnfun2(lambda,[0,1+dd,0,1]);
rhs = @(x,y) g(x,y);  % function g
[pts,al]=rational_kernel(order,'equi');
rhs2 = ultraSEM.Sol(rhs, n, dom);

%% Compute f
pdo={{@(x,y) -1+0*x+0*y,@(x,y) 0*x+0*y,@(x,y) -1+0*x+0*y},...
               {@(x,y) 0*x+0*y,@(x,y) 0*x+0*y},@(x,y) 0*x+0*y};
op=ultraSEM(dom, pdo, rhs, n); bc = 0;
f=op \ bc;

c = sqrt(sum2(f.*rhs2));
rhs = @(x,y) g(x,y)/c;  % normalise
rhs2 = (1/c)*rhs2;

%% Compute smoothed measure at evaluation points
mu=zeros(length(X),length(ep_vec));
pf = parfor_progress(length(X)*length(ep_vec)*order);
pfcleanup = onCleanup(@() delete(pf));

for lll=1:length(ep_vec)
	epsilon = ep_vec(lll);
    for i=1:length(X)
		I=0;
		for ii=1:order
			% solve systems using ultraSEM
			z=X(i)-pts(ii)*epsilon;
			pdo={{@(x,y) z+0*x+0*y,@(x,y) 0*x+0*y,@(x,y) z-1+0*x+0*y},...
				{@(x,y) 0*x+0*y,@(x,y) 0*x+0*y},@(x,y) 0*x+0*y};
			op=ultraSEM(dom, pdo, rhs, n); bc = 0;
			sol=op \ bc;        
			I=I+imag(al(ii)*sum2(sol.*rhs2));
			parfor_progress(pf);
		end
		mu(i,lll)=-I/pi;
		
    end
end
%%

figure
plot(X,mu)




function [poles,res] = rational_kernel(m,type)
%%rational_kernel 
% Inputs: 
%         - order of kernel 'm'
%         - Pole location 'roots','cheb'
% Outputs: - poles of rational function 'poles' 
%         - residues of rational function 'res'

%Poles
if type=="cheb"
    z=1i+chebpts(m,1);  %Chebyshev points
    %z=1i+(2*(1:m)/(m+1)-1);%.^2.*sign((2*(1:m)/(m+1)-1));
elseif type=="roots"
    z=exp(1i*pi*((1:m)'-1/2)/m); %Roots of unity
elseif type=="extrap"
    z=1i*(0.5.^(0:m-1));
elseif type=="equi"
    z=1i+(2*(1:m)/(m+1)-1);
end

if type=="equi" && m<7 %Hard-coded kernels from table of the paper
    if m==1
        res=1;
    elseif m==2
        res=[(1+3i)/2;(1-3i)/2];
    elseif m==3
        res=[-2+1i;5;-2-1i];
    elseif m==4
        res=[(-39-65i)/24;(17+85i)/8;(17-85i)/8;(-39+65i)/24];
    elseif m==5
        res=[(15-10i)/4;(-39+13i)/2;(65)/2;(-39-13i)/2;(15+10i)/4];
    else
        res=[(725+1015i)/(192);(-2775-6475i)/(192);(1073+7511i)/(96);(1073-7511i)/(96);(-2775+6475i)/(192);(725-1015i)/(192)];
    end
else %Vandermonde matrix for poles
    V=zeros(m);
    for i=1:m 
        V(:,i)=(z).^(i-1);
    end

    %Get residues directly
    rhs=eye(m,1);
    res=transpose(V)\rhs;
end

poles=z;

end

