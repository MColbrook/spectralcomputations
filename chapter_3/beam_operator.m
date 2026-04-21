clear
close all

%% Parameters

n1 = 100;
n2 = 50; % size of discretization is 2*n2+1
SCALE = 1; % scale factor

%% Set up the matrix

% Hermite basis
bb=spdiags(sqrt((0:n1+100)'/2),1,n1+100,n1+100);
mx=(bb+bb')/2; % multiplication by x/2
dx=bb-bb'; % derivative operator

mx2 = mx*mx; 
mx3 = mx*mx2; 
mx4 = mx*mx3; 

dx4 = dx*dx*dx*dx; dx4 = dx4(1:(n1+4),1:(n1+4));
mx = mx(1:(n1+4),1:(n1+4));
mx2 = mx2(1:(n1+4),1:(n1+4));
mx3 = mx3(1:(n1+4),1:(n1+4));
mx4 = mx4(1:(n1+4),1:(n1+4));


n2 = n2 + 10;
%  Malmquist-Takenaka basis
V1 = @(x) 1/(1+x.^2); % potential
V2 = @(x) x./(1+x.^2); % potential
cut_off = 10^(-30); % cut-off parameter to make operator banded

N = n2 + 30; % trick to get actual truncations of matrices
D1 = spdiags(transpose([-N:N+2; (-(N+1):(N+1))+1/2]),[-1,0],2*N+3,2*N+3);
D1 = 1i/2*(D1+D1')*SCALE;

D2 = D1*D1;
D2 = D2(2:end-1,2:end-1);
D2 = D2(31:end-30,31:end-30);
D3 = D1*D1*D1;
D3 = D3(2:end-1,2:end-1);
D3 = D3(31:end-30,31:end-30);
D4 = D1*D1*D1*D1;
D4 = D4(2:end-1,2:end-1);
D4 = D4(31:end-30,31:end-30);
D1 = D1(2:end-1,2:end-1);
D1 = D1(31:end-30,31:end-30);

v = chebfun(@(x) V1(1i*(1-exp(1i*x))./(SCALE*(1+exp(1i*x)))) ,[-pi,pi],'trig');
c = trigcoeffs(v,4*n2+1);
c(abs(c)<cut_off)=0;
V1 = sptoeplitz(c(2*n2+1:4*n2+1),c(2*n2+1:(-1):1));

v = chebfun(@(x) V2(1i*(1-exp(1i*x))./(SCALE*(1+exp(1i*x)))) ,[-pi,pi],'trig');
c = trigcoeffs(v,4*n2+1);
c(abs(c)<cut_off)=0;
V2 = sptoeplitz(c(2*n2+1:4*n2+1),c(2*n2+1:(-1):1));

n2 = n2 - 10;

%% Finite section matrix

E = [];
figure
for a = 0:0.1:3
    T = kron(dx4(1:n1,1:n1),speye(2*n2+1)) + kron(speye(n1),D4(11:end-10,11:end-10)) + 4i*kron(mx(1:n1,1:n1),D3(11:end-10,11:end-10))...
        - 6*kron(mx2(1:n1,1:n1),D2(11:end-10,11:end-10)) - 4i*kron(mx3(1:n1,1:n1),D1(11:end-10,11:end-10)) + kron(mx4(1:n1,1:n1),speye(2*n2+1))...
        + 2*a*kron(speye(n1),V2(11:end-10,11:end-10)) + a^2*kron(speye(n1),V1(11:end-10,11:end-10));
    E2 = real(eigs(T,100,'smallestabs'))';
    E = [E,E2+1i*a,E2-1i*a];
    plot(real(E),imag(E),'.k')
    xlim([-1,6])
    pause(0.1)
end

%% CompSpec - for speedup, change to computing residuals of finite sections

figure
for a = 0:0.1:3

    T = kron(dx4,speye(2*n2+21)) + kron(speye(n1+4),D4) + 4i*kron(mx,D3)...
            - 6*kron(mx2,D2) - 4i*kron(mx3,D1) + kron(mx4,speye(2*n2+21))...
            + 2*a*kron(speye(n1+4),V2) + a^2*kron(speye(n1+4),V1);
    
    % delete columns
    I = round(10000*log(kron(exp((1:n1+4)/10000),exp(1i*(-(n2+10):(n2+10))/10000))));
    J1 = find(real(I)>n1);
    J2 = find(abs(imag(I))>n2);
    
    Id = speye(size(T));
    T(:,[J1,J2]) = [];
    Id(:,[J1,J2]) = [];
    
    % faster way (can also do via residuals of finite section eigenpairs for even faster way)
    T2 = T; T2([J1,J2],:) = [];
    xpts = real(eigs(T2,100,'smallestabs'))';
    xpts = xpts(xpts<7);
    
    % slower way
    % xpts = -1:0.1:6;
    
    dist = 0*xpts;
    
    
    pf = parfor_progress(length(xpts));
    pfcleanup = onCleanup(@() delete(pf));
    for jj=1:length(xpts)
        B = T-xpts(jj)*Id;
        dist(jj) = sqrt(eigs(B'*B,1,'smallestabs'));
        parfor_progress(pf);
    end


    spec = 0*dist;
    for jj = 1:length(xpts)    
        if dist(jj)<=1
            x = find (abs(xpts-xpts(jj))<dist(jj));
            d=dist(x);
            spec(x(d==min(d(:))))=1;
        end
    end
 
    
    plot(xpts(spec==1),xpts(spec==1)*0+a,'.k')
    hold on
    plot(xpts(spec==1),xpts(spec==1)*0-a,'.k')
    xlim([-1,6])
    pause(0.1)
end






