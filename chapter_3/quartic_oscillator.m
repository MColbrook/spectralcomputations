clear
close all

% Add utils and data folders
currentFile = mfilename('fullpath');
[currentDir, ~, ~] = fileparts(currentFile);
targetPath = fullfile(currentDir, '..', 'utils');
addpath(genpath(targetPath));
targetPath = fullfile(currentDir, '..', 'data_online');
addpath(genpath(targetPath));

%%%%%%%%%%%%%%% CODE FOR QUARTIC OSCILLATOR EXAMPLE %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Set up the matrix
N=200;
bb=spdiags(sqrt((0:N+100)'/2),1,N+100,N+100);
mx=bb+bb';
dx=bb-bb';

dx=dx(1:(N+99),1:(N+99)); % derivative operator
mx=mx(1:(N+99),1:(N+99)); % multiplication by x

Lap=dx*dx;
T=-Lap+mx*mx*mx*mx;

EA = (3*gamma(3/4)*((1:N)-1/2)*sqrt(pi)/gamma(1/4)).^(4/3); % asymptotics of eigenvalues

%% Convergence of 5th eigenvalue

Nvec=10:100;
LL=length(Nvec);
ER = zeros(LL,1);
ER2=ER;
r = min(EA(6)-EA(5),EA(5)-EA(4))/4;

for jj=1:LL
    NN = Nvec(jj);
    [~,ER(jj)]=findmin(T(1:NN+4,1:NN),[EA(5)-r,EA(5)+r],10^(-13));

    B=T(1:NN+4,1:NN)-EA(5)*speye(NN+4,NN);
    [~,~,V]=svds(B,1,'smallest');

    rho = T(1:NN+4,1:NN)*V;
    rho = rho(1:length(V))'*V;
    epsilon = norm((T(1:NN+4,1:NN)-rho*speye(NN+4,NN))*V);
    r2 = min(abs(rho-11.7),abs(rho-21.1));

    ER2(jj) = epsilon^2/r2;
   
    figure(2)
    semilogy(Nvec,ER2,'.k','markersize',20) % in practice we estimate the denominator using distspec
    hold on
    semilogy(Nvec,ER,'ok','markersize',10)    
    hold off
    
    xlabel('$n$','interpreter','latex','fontsize',18)
    title('Comparison of a posteriori bounds','fontsize',18,'interpreter','latex')
    legend({'Bound on ${\|\!(A{-}\rho_{\mathcal{V}_n}\!I)x_{\mathcal{V}_n}\!\|^2}\!/{\delta(\rho_{\mathcal{V}_n}\!)}\!\,$','Local minimum of $\gamma_n$'},'fontsize',16,'interpreter','latex','location','southwest')
    
    xlim([10,100])
    ylim([10^(-15),10])
    ax = gca; ax.FontSize = 17;

    pause(0.001)

end


function [s,bd] = findmin(H2,s0,TOL)
Id2=speye(size(H2,1),size(H2,2));
while s0(2)-s0(1)>TOL   
    B=H2-mean(s0)*Id2;
    [~,~,V]=svds(B,1,'smallest');
    d1=norm(B*V);
    
    
    B=H2-mean(s0+TOL/2)*Id2;
    [~,~,V]=svds(B,1,'smallest');
    d2=norm(B*V);

    s=mean(s0);
    bd=d1;
    if d1-d2<0
        s0=[s0(1),mean(s0)];
    else
        s0=[mean(s0),s0(2)];
    end
end
end





