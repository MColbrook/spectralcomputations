function K = Kk (k,N)
%%
sm =  1;
for n=1:N
    pr = 1;
    for j=1:2:2*n-1
        pr=pr*j/(j+1);
    end
    sm = sm+pr^2*k^(2*n);
end
K = sm*pi/2;
%%
end