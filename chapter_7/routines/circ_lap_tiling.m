 function [L, type, loc, ang, indx_bdy] = circ_lap_tiling(level)

% function [L, type, loc, ang, indx_bdy] = cir_lap_tiling(level)
%
% Construct the Robinson stone version of the Penrose tiling
% and its graph Laplacian.
%
% Embree and Fillman, 2015-2022
% performance fine-tuned at MFO, August 2023
%
% This version starts from a 20-tile convex decagon seed.
% **
% ** type, loc, ang:  the type, location, and orientation of the tiles
% ** indx_bdy:        indices of the nodes on the boundary
% 

% types:  1 = B (blue)
%         2 = P (pink)
%         3 = W (white)
%         4 = Y (yellow)

% rule 1: B1 - P1:  
% rule 2: B2 - P2:
% rule 3: B3 - P3:
% rule 4: B2 - Y2:
% rule 5: B3 - W3:
% rule 6: P2 - W2:
% rule 7: P3 - Y3:
% rule 8: Y3 - W3:

% set the basic tile shapes
 phi = (1+sqrt(5))/2;

 for kk=0:level
    fprintf('working on level %d....\n', kk)

    if kk==0,   % initialize the first level
       types = kron([1;1;1;1;1],[1;3;2;4]);
       orig  = kron(exp(2i*pi*[1:5]'/5),[1;1;1;1]);

       edges = [ 1     2
                 3     4
                 5     6
                 7     8
                 9    10
                11    12
                13    14
                15    16
                17    18
                19    20
                 1     3
                 5     3
                 6     4
                 5     7
                 9     7
                10     8
                 9    11
                13    11
                14    12
                13    15
                17    15
                18    16
                17    19
                 1    19
                 2    20];

       rules = [ 5 7 5 7 5 7 5 7 5 7 3 2 8 3 2 8 3 2 8 3 2 8 3 2 8];
       th   = [2.5;1.3;2.5;2.7]*pi;
       angl  = mod(kron(2*pi*[1:5]'/5,[1;1;1;1])+kron([1;1;1;1;1],th),2*pi);
       N = length(types);
 
   else % compute the next level
        % all tiles at the previous level become 2 or 3 new tiles

      N = size(types,1);
      newtypes = NaN*ones(N,3);            % each new tile has a type (B,P,W,Y)
      newid    = NaN*ones(N,3);            % id of each new tile

      newedges = [];
      newrules = [];
      neworig = [];
      newangl = [];

% first let all tiles generate their spawn...
% and set the internal edges between the new tiles

      newedges = NaN*ones(2*N,2);
      id = 0;
     
      for j=1:N
         m = 3*(j-1)+1;   % refinement: below, no need for "newid" if you use "m"
         switch types(j)    % each B,P,W,Y always splits the same way 
            case 1, % B tiles -> P, Y, B tiles
                    newtypes(j,:) = [2 4 1];                     % B -> P Y B
                    newid(j,:)    = id+[1 2 3];  
                    newedges([2*j-1 2*j],[1 2]) = ...
                               [newid(j,1) newid(j,2)            % interior P <-> Y (7)
                                newid(j,3) newid(j,2)];          % interior B <-> Y (4)
                    newrules = [newrules; 7; 4]; 
                    neworig  = [neworig;phi*orig(j)+exp(1i*angl(j))*[1i;1i;phi*1i]];
                    newangl  = [newangl;angl(j)+pi*[1;6/5;4/5]];
                    id = id+3;
            case 2, % P tiles -> B, W, P tiles
                    newtypes(j,:) = [1 3 2];                     % P -> B W P
                    newid(j,:)    = id+[1 2 3];  
                    newedges([2*j-1 2*j],[1 2]) = ...
                               [newid(j,1) newid(j,2)            % interior B <-> W (5)
                                newid(j,3) newid(j,2)];          % interior P <-> W (6)
                    newrules = [newrules; 5; 6]; 
                    neworig  = [neworig;phi*orig(j)+exp(1i*angl(j))*[1i;1i;phi*1i]];
                    newangl  = [newangl;angl(j)+pi*[1;-1/5;-4/5]];
                    id = id+3;
            case 3, % W tiles -> W, B tiles
                    newtypes(j,:) = [3 1 NaN];                   % W -> W B
                    newid(j,:)    = id+[1 2 NaN];  
                    newedges(2*j-1,[1 2]) = ...
                               [newid(j,2) newid(j,1)];          % interior B <-> W (5)
                    newrules = [newrules; 5]; 
                    neworig  = [neworig;phi*orig(j)+exp(1i*angl(j))*exp(1i*11*pi/10)*[-1;-1]];
                    newangl  = [newangl;angl(j)+pi*[-3/5;3/5]];
                    id = id+2;
            case 4, % Y tiles -> Y, P tiles
                    newtypes(j,:) = [4 2 NaN];                   % Y -> Y P
                    newid(j,:)    = id+[1 2 NaN];  
                    newedges(2*j-1,[1 2]) = ...
                               [newid(j,2) newid(j,1)];          % interior P <-> Y (7)
                    newrules = [newrules; 7]; 
                    neworig  = [neworig;phi*orig(j)+exp(1i*angl(j))*exp(1i*9*pi/10)*[-1;-1]];
                    newangl  = [newangl;angl(j)+pi*[3/5;2/5]];
                    id = id+2;
         end
      end
      indx = find(~isnan(newedges(:,1)));
      newedges = newedges(indx,[1 2]);

% establish edges between tiles generated from different parents

% rule 1: B1 - P1: ->  2, 8
% rule 2: B2 - P2: ->  3
% rule 3: B3 - P3: ->  3, 1
% rule 4: B2 - Y2: ->  3
% rule 5: B1 - W2: ->  2, 8
% rule 6: P2 - W2: ->  3
% rule 7: P3 - Y3: ->  2, 8
% rule 8: W1 - Y1: ->  9
% rule 9: W3 - Y2: ->  3

      newpedges = NaN*ones(2*size(edges,1),2);
      for j=1:size(edges,1)
         j1 = edges(j,1);
         j2 = edges(j,2);
         switch rules(j)
             case 1, newpedges([2*j-1 2*j],[1 2]) = ...
                                [newid(j2,1) newid(j1,1)          % exterior B2 <-> P1 (2)
                                 newid(j2,2) newid(j1,2)];        % exterior W1 <-> Y2 (8)
                     newrules = [newrules; 2; 8];
             case 2, newpedges([2*j-1],[1 2]) = ...
                                [newid(j1,3) newid(j2,3)];        % exterior B1 <-> P2 (3)
                     newrules = [newrules; 3];
             case 3, newpedges([2*j-1 2*j],[1 2]) = ...
                                [newid(j2,1) newid(j1,1)          % exterior B2 <-> P1 (3)
                                 newid(j1,3) newid(j2,3)];        % exterior B1 <-> P2 (1)
                     newrules = [newrules; 3; 1];
             case 4, newpedges([2*j-1],[1 2]) = ...
                                [newid(j1,3) newid(j2,2)];        % exterior B1 <-> Y2 (3)
                     newrules = [newrules; 3];
             case 5, newpedges([2*j-1 2*j],[1 2]) = ...
                                [newid(j2,2) newid(j1,1)          % exterior B2 <-> P1 (2)
                                 newid(j2,1) newid(j1,2)];        % exterior W1 <-> Y2 (8)
                     newrules = [newrules; 2; 8];
             case 6, newpedges([2*j-1],[1 2]) = ...
                                [newid(j2,2) newid(j1,3)];        % exterior B2 <-> P1 (3)
                     newrules = [newrules; 3];
             case 7, newpedges([2*j-1 2*j],[1 2]) = ...
                                [newid(j1,1) newid(j2,2)          % exterior B1 <-> P2 (2)
                                 newid(j1,2) newid(j2,1)];        % exterior W1 <-> Y2 (8)
                     newrules = [newrules; 2; 8];
             case 8, newpedges([2*j-1],[1 2]) = ...
                                [newid(j1,1) newid(j2,1)];        % exterior W1 <-> Y2 (1)
                     newrules = [newrules; 9];
             case 9, newpedges([2*j-1],[1 2]) = ...
                                [newid(j1,2) newid(j2,2)];        % exterior B1 <-> P2 (3)
                     newrules = [newrules; 3];
         end
      end
      indx = find(~isnan(newpedges(:,1)));
      newedges = [newedges; newpedges(indx,[1 2])];
      types = reshape(newtypes',3*N,1);
      indx  = find(~isnan(types));
      types = types(indx);
      edges = newedges;
      rules = newrules;
      orig  = neworig;
      angl  = newangl;
   end
 end

 N = length(types);
 fprintf('Level %d:  %5d tiles\n', kk, N)

% create L; vast speed-up to create L with all the edges at once
% the other lines (L=L+L', L = diag(d)-L) are fast
% Note: replacing diag(d) with spdiags(d,0,N,N)
% significantly slows down the construction....

 edges = sort(edges,2);
 elen = size(edges,1);

 L = sparse(edges(:,1),edges(:,2),ones(elen,1),N,N,4*N);
 L = L+L';
 d = sum(L,2);     % d is a sparse vector...
 L =  diag(d)-L;   % and so diag(d) is a sparse matrix
 indx_bdy = find(d<3);

 type = types;
 loc  = orig;
 ang  = angl;
