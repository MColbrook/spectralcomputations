function [Yp,Yac] = mu_type(U,c,Nmvec,nvec,mvec)
% Computes spectral type (pure point part and absolutely continuous part) for the unitary operator U
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:       U: discretisation of unitary operator
%               c: coefficients of vector we compute spectral measure w.r.t.
%               Nmvec: the truncation sizes for Cesaro/Fourier sums
%               nvec: truncation parameters (n2) for pp part
%               mvec: truncation parameters (n2) for ac part
%                           
% OUTPUTS:      Yp: pure point part of total spectral measure
%               Yac: absolutely continuous part of total spectral measure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute point spectra part

N = (size(U,1)-1)/2;
Nm = max(Nmvec);
MU = zeros(Nm+1,1);     MU(1) = 1;

v = c;
MU1 = zeros(Nm,length(nvec));
for j=1:Nm
    v = U*v;
    MU(j+1) = c'*v;
    for kk=1:length(nvec)
        n = nvec(kk);
        II = (N+1-n):(N+1+n);
        MU1(j,kk) = norm(v(II));
    end
end

v = c;
MU2 = zeros(Nm,length(nvec));
for j=1:Nm
    v = U'*v;
    for kk=1:length(nvec)
        n = nvec(kk);
        II = (N+1-n):(N+1+n);
        MU2(j,kk) = norm(v(II));
    end
end

Yp = zeros(Nm,length(nvec));
for jj=Nmvec
    Yp(jj,:)=(1+sum(MU1(1:jj,:)+MU2(1:jj,:),1))/(2*jj+1);
end
Yp = Yp(Nmvec,:);

%% Compute ac spectra part

MU=flipud([conj(flipud(MU(2:end)));MU]/(2*pi));

Yac = zeros(length(Nmvec),length(mvec));
for jj=1:length(Nmvec)
    
    nn = abs(-Nmvec(jj):Nmvec(jj));
    gJ = ((1-nn/Nmvec(jj)).*cos(nn/Nmvec(jj)*pi) + cot(pi/Nmvec(jj))/Nmvec(jj)*sin(nn/Nmvec(jj)*pi));
    II = ((length(MU)-1)/2+1 - Nmvec(jj)):((length(MU)-1)/2+1 + Nmvec(jj));

    nu1 = chebfun(gJ(:).*MU(II),[-pi pi],'trig','coeffs');
    nu1 = real(nu1(-pi:0.0001:pi));
    
    for kk=1:length(mvec)
        II = nu1<mvec(kk);
        Yac(jj,kk) = sum(nu1(II))*0.0001;
    end

end



end