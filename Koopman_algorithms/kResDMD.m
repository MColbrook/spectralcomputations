function [G,K,L,PX,PY,PSI_x,PSI_y] = kResDMD(X,Y,varargin)
% This code applies kernelized ResDMD for the Perron--Frobenius operator on
% an RKHS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% Xa and Ya: data matrices used in kernel_EDMD to form dictionary (columns
% correspond to instances of the state variable)

% OPTIONAL LABELLED INPUTS
% N: size of computed dictionary, default is number of data points for kernel EDMD
% type: kernel used, default is normalised Gaussian, "Laplacian" is for
% normalised Laplacian, and numeric value (e.g., 20) is for polynomial
% kernel
% cut_off: stability parameter for SVD, default is 0

% OUTPUTS
% Ghat Khat Rhat matrices for kResDMD
% PSI matrices for evaluations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Collect the optional inputs
p = inputParser;

addParameter(p,'N',size(X,2),@(x) x==floor(x))
addParameter(p,'type',"Gaussian");
addParameter(p,'cut_off',10^(-14),@(x) x>=0)
addParameter(p,'kernel',[])
addParameter(p,'Xb',[],@isnumeric)
addParameter(p,'Yb',[],@isnumeric)

p.CaseSensitive = false;
parse(p,varargin{:})



%% Form the kernel

if ~isempty(p.Results.kernel)
    kernel_f = @(x,y) p.Results.kernel(x,y);
elseif isnumeric(p.Results.type)
    d = mean(vecnorm(X));
    kernel_f = @(x,y) transpose((y'*x/d^2+1).^(p.Results.type));
elseif p.Results.type=="Linear"
    kernel_f = @(x,y) transpose(y'*x);
elseif p.Results.type=="Laplacian"
    d = mean(vecnorm(X-mean(X,2)));
    if isa(X,'single') % safeguard against square root (but a little bit slower)
        kernel_f = @(x,y) transpose(exp(-pdist2(y',x')/d));
    else
        kernel_f = @(x,y) transpose(exp(-sqrt(-2*real(y'*x)+dot(x,x)+dot(y,y)')/d));
    end
elseif p.Results.type=="Gaussian"
    d = mean(vecnorm(X-mean(X,2)));
    kernel_f = @(x,y) transpose(exp(-(-2*real(y'*x)+dot(x,x)+dot(y,y)')/d^2));
elseif p.Results.type=="Lorentzian"
    d = mean(vecnorm(X-mean(X,2)));
    kernel_f = @(x,y) transpose((1+(-2*real(y'*x)+dot(x,x)+dot(y,y)')/d^2).^(-1));
end

%% Apply kernel EDMD

G1 = kernel_f(X,X); G1 = (G1+G1')/2;
A1 = kernel_f(Y,X);
L1 = kernel_f(Y,Y);  L1 = (L1+L1')/2;

% Post processing

[U,D0] = eig(G1+norm(G1)*p.Results.cut_off*eye(size(G1)));
[~,I] = sort(diag(D0),'descend');
U = U(:,I); D0 = D0(I,I);
N = min(p.Results.N,length(find(diag(D0)>0)));
U = U(:,1:N); D0 = D0(1:N,1:N);
UU = U*sqrt(diag(1./diag(D0)));

G = eye(N);
K = transpose(UU'*A1*UU);
L = transpose(UU'*L1*UU);

PX = G1*UU;
PY = A1*UU;

%% Evaluate on test data if wanted

if ~isempty(p.Results.Xb) % test data case
    PSI_x = kernel_f(p.Results.Xb,X)*UU;
else
    PSI_x =[];
end

if ~isempty(p.Results.Yb) % test data case
    PSI_y = kernel_f(p.Results.Yb,X)*UU;
else
    PSI_y =[];
end

end
