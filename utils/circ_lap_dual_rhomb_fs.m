  function [L, loc, edges, indx_bdy] = circ_lap_dual_rhomb_fs(level)

% Generate graph Laplacian for *finite section of the infinite tiling*.
% NB:  Generates Laplacian for level+2 to determine the degree of the 
% boundary nodes in the original tiling.
%
% February 2024

 k = level;
 [L,loc,edges,indx_bdy] = circ_lap_dual_rhomb(k);
 [L2,loc2,edges2,indx_bdy2] = circ_lap_dual_rhomb(k+2);

 loc2 = -loc2;     % rotate 180 degrees, i.e., multiply by exp(pi*i)=-1

 N = length(loc);
 N2 = length(loc2);

 indx = NaN*ones(N,1);
 d    = NaN*ones(N,1);
 for mm=1:length(indx_bdy)
    m = indx_bdy(mm);
    [dist,j] = min(abs(loc2-loc(m)));
    if dist<1e-10, indx(m) = j; else, fprintf('ERROR: no good match for node %d\n',m); end 
    d(m) = L2(j,j);
    dold = full(L(m,m));
    L(m,m) = d(m);
    fprintf(' %2d;  old degree = %2d;  new degree = %2d\n', m, dold, d(m));
 end
