clear
close all

B = 1; % strength of magnetic field

for N = 1:50
    % Set up the matrix using recurrence relations for Hermite functions
    bb = spdiags(sqrt((0:N+100)'/(2)),1,N+100,N+100);
    mx = bb+bb';
    mx2 = mx*mx;
    dx = bb-bb';
    dx2 = dx*dx;
    
    mx = mx(1:N,1:N);
    mx2 = mx2(1:N,1:N);
    dx = dx(1:N,1:N);
    dx2 = dx2(1:N,1:N);
    
    I = speye(N);
    
    H = kron(-dx2,I)+1i*kron(dx,mx)*B+B^2/4*kron(I,mx2)...
        +kron(I,-dx2)-1i*kron(mx,dx)*B+B^2/4*kron(mx2,I);
    
    % Compute and plot eigenvalues
    E = eig(full(H));
    
    plot(N+0*E,E,'k.','markersize',12)
    hold on
    ylim([0,10])
    pause(0.01)
    if N==1 % plot the actual eigenvalues
        plot([0,1]*3,[1,1],'g','linewidth',3)
        plot([0,1]*3,[1,1]*3,'g','linewidth',3)
        plot([0,1]*3,[1,1]*5,'g','linewidth',3)
        plot([0,1]*3,[1,1]*7,'g','linewidth',3)
        plot([0,1]*3,[1,1]*9,'g','linewidth',3)
    end
end

% label things
xlabel('$N$ (discretisation size is $N^2$)','interpreter','latex','fontsize',18)
ylabel('Computed Eigenvalues','interpreter','latex','fontsize',18)
legend({'Finite Section Eigenvalues','True Eigenvalues'},'fontsize',16,'interpreter','latex','location','northwest')
ax=gca; ax.FontSize=18;
exportgraphics(gcf,'CHAP1_magpoll_hermite.pdf','ContentType','vector','BackgroundColor','none','Resolution',300,'Colorspace','gray')
