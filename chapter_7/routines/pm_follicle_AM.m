 function fig = pm_follicle_AM(sigma, q, plot_bdy, lw, col)

% function fig = pm_follicle(comp_spec, lambda_min, lambda_max, num_lambda, plot_bdy, lw, col)
%
% "follicle plot" showing the evolution of the spectrum as the 
% coupling constant lambda increases.
%
% fig        = handle to the figure 
% comp_spec  = handle to a function that computes the spectrum as real intervals;
%              - assumes one argument, the coupling constant lambda
%              - returns the spectrum as a K-by-2 matrix to describe K intervalsha
% lambda_min = minimum coupling constant (default = 0)
% lambda_max = maximum coupling constant (default = 2)
% num_lambda = number of lambda values (default = 1000)
% plot_bdy   = plot boundary as black dots at the end of each interval? 
%              0 => no boundary; > 0 => markersize of boundary (default = 0.5)
% lw         = linewidth (default = 450/num_lambda)
% col        = color (default [.7 .9 1])
%
% When printing to .pdf, it is best to use the "-painters" flag to get
% vector graphics output, e.g.
%    print -dpdf -painters foo1.pdf


% For each lambda value, the spectrum is drawn as a union of horizontal 
% intervals at height lambda.  The idea is to use sufficiently many lambda
% values that the resulting plot looks like a filled-in area, rather than 
% discrete lines.  This depends on the number of lambda values and width of
% the lines used to draw the intervals.

 if nargin<3, plot_bdy = 0.1; end
 if nargin<4, lw = 255/q; end
 if nargin<5, col = [1,1,1]*0.9; end
     %[.7 .9 1]; end

 % lamvec = flipud(linspace(lambda_min,lambda_max,num_lambda)');
 fig = figure(1); clf
 minspec = Inf;
 maxspec = -Inf;

 dlam = 1/(2*q);
 for j=0:q
    S = sigma{j+1};%comp_spec(lamvec(j));                    % compute the itervals
    minspec = min([min(S(:,1)) minspec]);
    maxspec = max([max(S(:,2)) maxspec]);
    if iscell(S), S = S{end}; end                % extract last component in cell array
    Sp = reshape([S NaN*ones(size(S,1),1)]',1,3*size(S,1))+j/q*1i;
    if j>0
        plot(real(Sp),imag(Sp),'-','linewidth',lw,'color',col)
    else
        plot(real(Sp),imag(Sp)+1/(4*q),'-','linewidth',lw/2,'color',col)
        
    end

    
    hold on
    if plot_bdy>0
       Sr = reshape(S,prod(size(S)),1);
       Sr = kron(Sr,[1,1,NaN]')+1i*kron(ones(length(Sr),1),(j/q+2*(dlam/2)*[-1 1 NaN]'));
       plot(real(Sr),imag(Sr),'k-','linewidth',.1)
    end
 end

 ylim([0 1])
 dspec = (maxspec-minspec)*.075;
 xlim([minspec-dspec,maxspec+dspec])
 ax = gca;
ax.Layer = 'top';
