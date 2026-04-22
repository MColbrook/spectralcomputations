function dim = box_dim(sigma,delta,er)
%%%%%% INPUTS
% sigma: a finite collection of intervals. First column is left endpoints, second column is right endpoints.
% delta: vector of covering lengths.
% er: pad the intervals with error tolerances if needed (default scales as smallest delta if unspecified)

%%%%%% OUTPUT
% dim: matrix of approximations of dim_B

if nargin<3
    er=min(delta)/10;
    tol = 10^(-16);
end

[~,I] = sort(sigma(:,1),'ascend');
sigma = sigma(I,:);

%% Process sigma

if er>0
    sigma(:,1) = sigma(:,1) - er; % pad with error if needed
    sigma(:,2) = sigma(:,2) + er;
    sigma=consolidate_int(sigma,tol);
end

%% Compute number of covers

dim = zeros(size(delta));
ee = rand(1,1);
sigma = sigma + ee;
for jj=1:length(delta)
    sig(:,1)=floor(sigma(:,1)/delta(jj));
    sig(:,2)=ceil(sigma(:,2)/delta(jj));

    [~,I]=sort(sig(:,1),'ascend');
    sig=sig(I,:);

    I = sig(1,2)-sig(1,1) + 1;
    m = sig(1,2);

    for j=2:size(sig,1)
        if sig(j,2)>m
            if sig(j,1)>m
                I = I + sig(j,2)-sig(j,1) + 1;
            else
                I = I + sig(j,2)-m;
            end
            m = sig(j,2);
        end
    end

    dim(jj) = log(I)/log(1/delta(jj));

end

end
