clear
close all

%%%%%% CODE FOR RESONANCE CONTOUR EXAMPLE (SECTION 10.3.5) %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters (some physical parameters hardcoded into nlevp_solver)

potential = @(x) 1+1/4-x/4+0.5*sin(3*(x+1)*2*pi); % continuous but constant for sufficiently large |x| (in this case \geq 1), this is it's value on [-1,1]
ksq = 121;
N = 300;                % number of Legendre polynomials for discretisation
nquad = 12:12:408;      % number of quadrature points for convergence plot
l = 75;                 % number of left probes
rr = 75;                % number of right probes
p = 1;                  % number of moments used
mreal = 64;
wght = @(bb) 1./(abs(bb)); % decay for the Gaussian process, needed in infinite dimensions (see spectral pollution example)

%% Gauss quadrature on each side of rectangle

% Random probing matrices
rng(0);
F = diag(wght(1:N))*randn(N,l);
G = diag(wght(1:N))*randn(N,rr);

delta = 0.1;
ht = 10;
Z0 = mean([delta,0.6])+mean([-0.1,ht])*1i;

ct = 1;
for n = nquad
    n
    [s1,w1] = legpts(round(n/4),[delta,0.6]);
    s1 = s1-0.1i;
    [s2,w2] = legpts(round(n/4),[-0.1,ht]);
    w2 = 1i*w2; s2 = s2*1i+0.6;
    [s3,w3] = legpts(round(n/4),[-0.1,ht]);
    w3 = -1i*w3; s3 = s3*1i+delta;
    [s4,w4] = legpts(round(n/4),[delta,0.6]);
    w4 = -w4; s4 = s4+ht*1i;
    
    poles = [s1(:);s2(:);s3(:);s4(:)];
    residues = [w1(:);w2(:);w3(:);w4(:)]/(2*pi*1i);
    
    A = cell({2*p});
    for jj = 0:2*p-1
        A{jj+1} = zeros(l,rr);
    end
    
    n = length(poles(:));
    for j = 1:n
        U = F'*nlevp_solver(poles(j), G, ksq, potential);
        for jj = 0:2*p-1
            A{jj+1} = A{jj+1} + residues(j)*(poles(j)-Z0)^jj*U/2i/pi;
        end
    end
    
    B0 = zeros(p*l,p*rr); B1 = B0;
    for aa = 1:p
        for bb = 1:p
            B0((l*(aa-1)+1):aa*l,(rr*(bb-1)+1):bb*rr) = A{aa+bb-1};
            B1((l*(aa-1)+1):aa*l,(rr*(bb-1)+1):bb*rr) = A{aa+bb};
        end
    end
    
    [V0, S0, W0] = svd(B0, 0);
    V0 = V0(:,1:mreal);  S0 = S0(1:mreal,1:mreal); W0 = W0(:,1:mreal);
    [~, ee] = eig(V0'*B1*W0,S0,'vector');
    [~,I] = sort(imag(ee)); ee = ee(I);
    lam_gauss{ct} = ee + Z0;
    ct = ct + 1;

end


%% Run trapezoidal rule individually around each eigenvalue

l = 2; % number of left probes
rr = 2; % number of right probes

nquad2=1:15;
clear lam
LAM = lam_gauss{end};
for m = 1:mreal
    ct = 1;
    m
    F = diag(wght(1:N))*randn(N,l);
    G = diag(wght(1:N))*randn(N,rr);
    for n = nquad2
        
        % small contour around each individual eigenvalue (much faster and
        % more accurate)
        cntr = LAM(m);
        Z0 = cntr;
        LL = 0.01;
        contour = @(t) cntr + LL*cos(t)+LL*1i*sin(t);
        jacobian = @(t) -LL*sin(t)+LL*1i*cos(t);

        tpts = linspace(0,2*pi,n+1)+2*pi/(2*n); tpts(end)=[];  % trap rule
        poles = contour(tpts);
        quadwts = (2*pi)/n * ones(1,n);
        residues = quadwts(:).*jacobian(tpts(:));



        A = cell({2*p});
        for jj = 0:2*p-1
            A{jj+1} = zeros(l,rr);
        end
        
        n = length(poles(:));
        for j = 1:n
            U = F'*nlevp_solver(poles(j), G, ksq, potential);
            for jj = 0:2*p-1
                A{jj+1} = A{jj+1} + residues(j)*(poles(j)-Z0)^jj*U/2i/pi;
            end
        end
        
        B0 = zeros(p*l,p*rr); B1 = B0;
        for aa = 1:p
            for bb = 1:p
                B0((l*(aa-1)+1):aa*l,(rr*(bb-1)+1):bb*rr) = A{aa+bb-1};
                B1((l*(aa-1)+1):aa*l,(rr*(bb-1)+1):bb*rr) = A{aa+bb};
            end
        end
        
        [V0, S0, W0] = svd(B0, 0);
        V0 = V0(:,1);  S0 = S0(1,1); W0 = W0(:,1);
        [~, ee] = eig(V0'*B1*W0,S0,'vector');
        [~,I] = sort(imag(ee)); ee = ee(I);
        lam(m,ct) = ee + Z0;
        ct = ct+1;
    end
end

%%

for ct = 1:length(nquad)-1
    E1(ct) = max(abs(lam_gauss{ct}-lam(:,end)));
end
%% Plot the results

figure
semilogy(nquad(1:end-1)*75,E1,'k.-','linewidth',2,'markersize',22)
xlim([0,30000])
ylim([10^(-15),10])
grid on
xlabel('total number of system solves','interpreter','latex','fontsize',18)
ylabel('max eigenvalue error','interpreter','latex','fontsize',18)
title('Single-Contour Gauss--Legendre','interpreter','latex','fontsize',18)
ax =gca;
ax.XAxis.Exponent = 0;
ax=gca; ax.FontSize=17;


figure
semilogy(nquad2*2*length(LAM),max(abs(lam(:,end)-lam),[],1),'k.-','linewidth',2,'markersize',22)
xlim([0,2000])
ylim([10^(-15),10])
grid on
xlabel('total number of system solves','interpreter','latex','fontsize',18)
ylabel('max eigenvalue error','interpreter','latex','fontsize',18)
title('Multiple-Contour Trapezoidal','interpreter','latex','fontsize',18)
ax =gca;
ax.XAxis.Exponent = 0;
ax=gca; ax.FontSize=17;


%% Solvers for T(z)^{-1}

function Sol = nlevp_solver(z, RHS, ksq, potential)
RHS = leg2cheb(RHS, 'norm');
N = size(RHS,1);
D2 = ultraS.diffmat(N, 2); 
S02 = ultraS.convertmat(N, 0, 1);

nl = potential(-1); nr = potential(1);

% careful branch cuts
ang = -1/4;
alphal = @(z) -sqrt(ksq)*sqrt(exp(2*pi*1i*ang)*(nl-z))*exp(-2*pi*1i*ang/2).*sqrt(nl+z);
alphar = @(z) -sqrt(ksq)*sqrt(exp(2*pi*1i*ang)*(nr-z))*exp(-2*pi*1i*ang/2).*sqrt(nr+z);

eta_sq = ultraS.multmat(N,chebfun(@(x) potential(x)),2);
T = D2 + ksq*(eta_sq*S02-z.^2*S02);

B=[1.^(0:N-1); (-1).^(0:N-1)];
for ii=1:1
    b=(0:N-1).^2;
    for j=1:ii-1
        b=b.*(((0:N-1).^2-j^2)/(2*j+1));
    end
    B=[B; b; (-1).^(ii+(0:N-1)).*b];
end


BC = [B(4,:)-1i*alphal(z)*B(2,:);
      B(3,:)+1i*alphar(z)*B(1,:)];

T = [speye(2,N);T(1:end-2,:)];
            
RHSb = S02*RHS;
RHSb = [zeros(2,size(RHSb,2)); RHSb(1:(N-2),:)];

U = eye(N,2); V = BC-eye(2,N);

Q1 = T\RHSb;
Q2 = T\U;
Sol = Q1-Q2*((eye(2,2)+V*Q2)\(V*Q1));
Sol = cheb2leg(Sol, 'norm');
end


