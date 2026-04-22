function dim = Haus_dim(sigma,n1vec,n2vec,er,cut)

%%%%%% INPUTS
% sigma: a finite collection of intervals. First column is left endpoints, second column is right endpoints.
% n1vec: n1 values for the tower of algorithms
% n2vec: n2 values for the tower of algorithms
% er: pad the intervals with error tolerances if needed (default is zero if unspecified)
% cut: point at which to measure increase in H^d (default is 1 if unspecified)

%%%%%% OUTPUT
% dim: matrix of approximations of dim_H

if nargin<4
    er = 0;
    cut = 1;
elseif nargin<5
    cut = 1;
end

%% Process sigma
sigma(:,1) = sigma(:,1) - er; % pad with error if needed
sigma(:,2) = sigma(:,2) + er;
sigma = consolidate_int(sigma,er);

sigma = sigma - min(sigma(:));
sigma_max = ceil(max(sigma(:)));

%% Compute h functions
dim = zeros(length(n1vec),length(n2vec));
for j=1:length(n1vec)
    n1 = n1vec(j) % print this out to keep track

    sig = zeros(size(sigma));
    sig(:,1) = floor(sigma(:,1)*2^(n1));
    sig(:,2) = ceil(sigma(:,2)*2^(n1));
    
    a = sig(:,1).'+1;
    b = sig(:,2).';
    cover = false(1,sigma_max*2^(n1));

    % Create an array of indices to be set to true
    indices = arrayfun(@(start, stop) start:stop, a, b, 'UniformOutput', false);

    % Concatenate the indices arrays
    indices = [indices{:}];

    cover(indices) = true;


    for kk=1:length(n2vec)
        
        n2 = n2vec(kk);
        n3 = n2;
        if n2vec(kk)>n1
            dim(j,kk)=NaN;
        else
            h=min_cover(cover,n1,n2,1);
            if h>=cut
                dim(j,kk)=1;
            else
                d = 1/2;
                ct = 1;
                while ct<n3+1
                    h=min_cover(cover,n1,n2,d);
                    if h<cut
                        d = d - 2^(-ct);
                    else
                        d = d + 2^(-ct);
                    end
                    ct = ct + 1;
                end
                h=min_cover(cover,n1,n2,d);
                while (h<cut)&&(d>0)
                    d=d-2^(-n3);
                    h=min_cover(cover,n1,n2,d);
                end
				max(min(d,1),0)
                dim(j,kk)=max(min(d,1),0);
            end
        end
    end
    figure(2)
    plot(n1vec,dim,'.-','linewidth',1,'markersize',12)
    pause(0.1)    
end

end



function [ h ] = min_cover(cover,n1,n2,d)

nsplit = round(length(cover)*2^(n2-n1));
cover = (reshape(cover,2^(n1-n2),nsplit));
cover = cover(:,any(cover,1));
h = zeros(1,size(cover,2));

for j = 1:size(cover,2)
    newcover = cover(:,j);
    h(j) = levelDP(n1-n2,d,n2,newcover);
end
h = sum(h);
end

function p = levelDP(ll,d,n,c)
    if ll == 0
        if any(c)
            p=2^(-d*n);
        else
            p = 0;
        end
    else
        b = 0;
        if any(c(1:2^(ll-1)))
            b = b + levelDP(ll-1,d,n+1,c(1:2^(ll-1)));
        end
        if any(c((2^(ll-1)+1):2^ll))
            b = b + levelDP(ll-1,d,n+1,c((2^(ll-1)+1):2^ll));
        end
        p=min(2^(-n*d),b);
    end
end

