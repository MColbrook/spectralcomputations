function [H,f,E] = ham_construct(phi,NN)
%H_HFST Hamiltonain for HS-model

% Hofstader model Hamiltonian
% This is a rather slow code for constructing the Hamiltonian! There are
% faster methods for large NN.

J = 1;                              %tunnelling strength
A_f = 1.023*0.7434*0.5;             %area factor
load('lattice_5Fold_117998.mat');   %load lattice structure
    
ltConnect(NN+1:end)=[];
ltPos_indX(NN+1:end)=[];
ltPos_indY(NN+1:end)=[];
f=zeros(NN,1);
E=ltPos_indX+1i*ltPos_indY;

for ii=1:length(ltPos_indX)
    ltConnect{ii}=ltConnect{ii}(ltConnect{ii}<NN+1);
    NNSitesMap(ii) = length(ltConnect{ii});
    for iii=1:length(ltConnect{ii})
        J_List{ii}(iii) = J;
    end
    f(ii)=max(ltConnect{ii});
end

nsites = length(ltConnect);
H = sparse(nsites,nsites);

H = H_hfst_SP(H,{J_List,A_f,phi},ltConnect,ltPos_indX,ltPos_indY);         % single particle hamiltonian


end

function H = H_hfst_SP(H,constList,ltConnect,ltX,ltY)
%H_HFST Hamiltonain for HS-model

J_List = constList{1};
A_f = constList{2};
phi = constList{3};

for i=1:length(ltConnect)
    for j=1:length(ltConnect{i})
        x_i = ltX(i);   y_i = ltY(i);
        x_j = ltX(ltConnect{i}(j)); y_j = ltY(ltConnect{i}(j)); 
        theta_ij = 2*pi*phi*(x_i + x_j)*(y_j-y_i)/(2*A_f);
        H(i,ltConnect{i}(j)) = -J_List{i}(j)*exp(1i*theta_ij);
    end
end

end


