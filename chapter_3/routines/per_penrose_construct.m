function [A] = per_penrose_construct(n)

%% constructs periodic approximations using the multigrid method
% NB: results may differ from run to run due to a random parameter used in
% the construction. The code is also not optimised.

% set the parameters
gamma = zeros(5,1);
gamma(2:4) = 0.1*rand(3,1); % random parameter ensure the grid is regular almost surely
gamma(5) = -sum(gamma(1:4));

small_param=0.00001;
F = @(n) (((1+sqrt(5))/2)^(n+1)-(-(1+sqrt(5))/2)^(-(n+1)))/sqrt(5);
tau=F(n+1)/F(n);


vects=zeros(5,2);
vects_perp=zeros(5,2);
phi=2*pi/5;
kvec=-2*F(n+1):2*F(n+1);

for j=1:5
    vects(j,:)=[cos((j-1)*phi),sin((j-1)*phi)];
    vects_perp(j,:)=[-sin((j-1)*phi),cos((j-1)*phi)];
end

%generate the approximants
vects_apr=vects;
vects_apr(3,:)=-vects(1,:)+(1/tau)*vects(2,:);
vects_apr(4,:)=-(1/tau)*(vects(1,:)+vects(2,:));
vects_apr(5,:)=(1/tau)*vects(1,:)-vects(2,:);

for j=1:5

    vects_perp(j,:)=[-vects_apr(j,2),vects_apr(j,1)];
end


% figure(1)
% hold on
p=F(n+1)*vects(2,:)+F(n)*vects(3,:)-F(n)*vects(4,:)-F(n+1)*vects(5,:);
q=-F(n+1)*vects(1,:)+F(n+1)*vects(3,:)+F(n)*vects(4,:)-F(n)*vects(5,:);
c=(p+q)/2;
vertices=[];


% plot([0, p(1),p(1)+q(1),q(1),0]-c(1),[0, p(2),p(2)+q(2),q(2),0]-c(2),'r')
% plot([p(1),q(1)]-c(1),[p(2),q(2)]-c(2),'ok','markersize',10)

% Now plot corner points to check code is working
for j=1:5
    z2=vects_apr(j,:)*vects_apr(j,:)';

    for l1=1:length(kvec)
        for k=1:5
            if k~=j
                z1=vects_apr(j,:)*vects_apr(k,:)';
                
                z3=vects_apr(k,:)*vects_apr(k,:)';
                for l2=1:length(kvec)
                    alpha=((kvec(l1)+gamma(j))*z3-z1*(kvec(l2)+gamma(k)))/(z2*z3-z1^2);
                    beta=(kvec(l1)+gamma(j)-alpha*z2)/(z1);
                    x=alpha*vects_apr(j,:)+beta*vects_apr(k,:);
                    
                    % Now add the other points
                    x1=x+small_param*(vects_perp(j,:)+vects_perp(k,:))/2;
                    x2=x+small_param*(vects_perp(j,:)-vects_perp(k,:))/2;
                    x3=x+small_param*(-vects_perp(j,:)-vects_perp(k,:))/2;
                    x4=x+small_param*(-vects_perp(j,:)+vects_perp(k,:))/2;

                    kappa=zeros(5,4);
                    for g=1:5
                        kappa(g,1)=floor(x1*vects_apr(g,:)'-gamma(g));
                        kappa(g,2)=floor(x2*vects_apr(g,:)'-gamma(g));
                        kappa(g,3)=floor(x3*vects_apr(g,:)'-gamma(g));
                        kappa(g,4)=floor(x4*vects_apr(g,:)'-gamma(g));
                    end
%                     
                    y1=kappa(1,1)*vects(1,:)+kappa(2,1)*vects(2,:)+kappa(3,1)*vects(3,:)+kappa(4,1)*vects(4,:)+kappa(5,1)*vects(5,:);
                    y2=kappa(1,2)*vects(1,:)+kappa(2,2)*vects(2,:)+kappa(3,2)*vects(3,:)+kappa(4,2)*vects(4,:)+kappa(5,2)*vects(5,:);
                    y3=kappa(1,3)*vects(1,:)+kappa(2,3)*vects(2,:)+kappa(3,3)*vects(3,:)+kappa(4,3)*vects(4,:)+kappa(5,3)*vects(5,:);
                    y4=kappa(1,4)*vects(1,:)+kappa(2,4)*vects(2,:)+kappa(3,4)*vects(3,:)+kappa(4,4)*vects(4,:)+kappa(5,4)*vects(5,:);
                    
                    center=(y1+y2+y3+y4)/4;
                    % Now add warning system if small_param is too large
                    check=0;
                    for g=1:5
                        if g~=j
                            if g~=k
                                check=check + (kappa(g,1)-kappa(g,2))^2+(kappa(g,1)-kappa(g,3))^2+(kappa(g,1)-kappa(g,4))^2;
                            end
                        end
                    end
                    if check>0
                        fprintf('ERROR: small_param too large');
                        return
                    end
                    I=inpolygon(center(1),center(2),[0, p(1),p(1)+q(1),q(1),0]-c(1),[0, p(2),p(2)+q(2),q(2),0]-c(2));
                    if I==1
                        vertices=[vertices;y1;y2;y3;y4];
%                         plot([y1(1),y2(1),y3(1),y4(1),y1(1)],[y1(2),y2(2),y3(2),y4(2),y1(2)],'k')

                    end
                end
            end
        end
    end
end
                 
% axis equal


% Now get rid of repetitions of vertices
I=size(vertices);
L=I(1);
I=[];

for a=1:L
    for b=a+1:L
        if (sum((vertices(a,:)-vertices(b,:)).^2)<0.0000001)
            I=[I;b];
        end
    end
end


V=vertices;
V(I,:)=[];
I2=size(V);
L=I2(1);

% for j=1:L
%     plot(V(j,1),V(j,2),'xb')
% end

% Now find the boundary vertices
J1=[];
J2=[];
J3=[];
J4=[];
for j=1:L
    if min((V(:,1)-V(j,1)+p(1)).^2+(V(:,2)-V(j,2)+p(2)).^2)<0.000001
        J1=[J1;j];
%         plot(V(j,1),V(j,2),'og')
    end
    if min((V(:,1)-V(j,1)+q(1)).^2+(V(:,2)-V(j,2)+q(2)).^2)<0.000001
        J2=[J2;j];
%         plot(V(j,1),V(j,2),'om')
    end
    if min((V(:,1)-V(j,1)+q(1)-p(1)).^2+(V(:,2)-V(j,2)+q(2)-p(2)).^2)<0.000001
        J3=[J3;j];
%         plot(V(j,1),V(j,2),'oy')
    end
    if min((V(:,1)-V(j,1)+q(1)+p(1)).^2+(V(:,2)-V(j,2)+q(2)+p(2)).^2)<0.000001
        J4=[J4;j];
%         plot(V(j,1),V(j,2),'oy')
    end
    
    
end

Vcell=V;
Vcell([J1; J2;J3;J4],:)=[]; % this is the set of (non repeating) vertices in our periodic cell

% now compute the connectivity matrix for periodic BCs

I2=size(Vcell);
L=I2(1);

Con_MAT_per=sparse(L,L);
for i=1:L
    I3=find((vertices(:,1)-Vcell(i,1)).^2+(vertices(:,2)-Vcell(i,2)).^2<0.0000001);
    L2=length(I3);
    for j=1:L2
        if mod(I3(j),4)==0
            v1=I3(j)-1; % the vs correspond to indices we join in the connectivity matrix
            
            v2=I3(j)-3;
            
        elseif mod(I3(j),4)==1
            v1=I3(j)+1;
            
            v2=I3(j)+3;
            
        elseif mod(I3(j),4)==2
            v1=I3(j)-1;
            v2=I3(j)+1;
            
        else
            
            v1=I3(j)-1;
            v2=I3(j)+1;
        end
        find1=find((vertices(v1,1)-Vcell(:,1)).^2+(vertices(v1,2)-Vcell(:,2)).^2<0.0000001);
        find2=find((vertices(v2,1)-Vcell(:,1)).^2+(vertices(v2,2)-Vcell(:,2)).^2<0.0000001);
        
        if isempty(find1)==0
            Con_MAT_per(find1,i)=1;
            Con_MAT_per(i,find1)=1;
        end
        if isempty(find2)==0
            Con_MAT_per(find2,i)=1;
            Con_MAT_per(i,find2)=1;
        end
        
       
    end
end

% Now we need to include boundary points
for i=1:L
    I3=find((vertices(:,1)-Vcell(i,1)-p(1)).^2+(vertices(:,2)-Vcell(i,2)-p(2)).^2<0.0000001);
    L2=length(I3);
    for j=1:L2
        if mod(I3(j),4)==0
            v1=I3(j)-1; % the vs correspond to indices we join in the connectivity matrix
            
            v2=I3(j)-3;
            
        elseif mod(I3(j),4)==1
            v1=I3(j)+1;
            
            v2=I3(j)+3;
            
        elseif mod(I3(j),4)==2
            v1=I3(j)-1;
            v2=I3(j)+1;
            
        else
            
            v1=I3(j)-1;
            v2=I3(j)+1;
        end
        find1=find((vertices(v1,1)-Vcell(:,1)).^2+(vertices(v1,2)-Vcell(:,2)).^2<0.0000001);
        find2=find((vertices(v2,1)-Vcell(:,1)).^2+(vertices(v2,2)-Vcell(:,2)).^2<0.0000001);
        
        if isempty(find1)==0
            Con_MAT_per(find1,i)=1;
            Con_MAT_per(i,find1)=1;
        end
        if isempty(find2)==0
            Con_MAT_per(find2,i)=1;
            Con_MAT_per(i,find2)=1;
        end
        
       
    end
end

for i=1:L
    I3=find((vertices(:,1)-Vcell(i,1)-q(1)).^2+(vertices(:,2)-Vcell(i,2)-q(2)).^2<0.0000001);
    L2=length(I3);
    for j=1:L2
        if mod(I3(j),4)==0
            v1=I3(j)-1; % the vs correspond to indices we join in the connectivity matrix
            
            v2=I3(j)-3;
            
        elseif mod(I3(j),4)==1
            v1=I3(j)+1;
            
            v2=I3(j)+3;
            
        elseif mod(I3(j),4)==2
            v1=I3(j)-1;
            v2=I3(j)+1;
            
        else
            
            v1=I3(j)-1;
            v2=I3(j)+1;
        end
        find1=find((vertices(v1,1)-Vcell(:,1)).^2+(vertices(v1,2)-Vcell(:,2)).^2<0.0000001);
        find2=find((vertices(v2,1)-Vcell(:,1)).^2+(vertices(v2,2)-Vcell(:,2)).^2<0.0000001);
        
        if isempty(find1)==0
            Con_MAT_per(find1,i)=1;
            Con_MAT_per(i,find1)=1;
        end
        if isempty(find2)==0
            Con_MAT_per(find2,i)=1;
            Con_MAT_per(i,find2)=1;
        end
        
       
    end
end



for i=1:L
    I3=find((vertices(:,1)-Vcell(i,1)-(p(1)+q(1))).^2+(vertices(:,2)-Vcell(i,2)-(p(2)+q(2))).^2<0.0000001);
    L2=length(I3);
    for j=1:L2
        if mod(I3(j),4)==0
            v1=I3(j)-1; % the vs correspond to indices we join in the connectivity matrix
            
            v2=I3(j)-3;
            
        elseif mod(I3(j),4)==1
            v1=I3(j)+1;
            
            v2=I3(j)+3;
            
        elseif mod(I3(j),4)==2
            v1=I3(j)-1;
            v2=I3(j)+1;
            
        else
            
            v1=I3(j)-1;
            v2=I3(j)+1;
        end
        find1=find((vertices(v1,1)-Vcell(:,1)).^2+(vertices(v1,2)-Vcell(:,2)).^2<0.0000001);
        find2=find((vertices(v2,1)-Vcell(:,1)).^2+(vertices(v2,2)-Vcell(:,2)).^2<0.0000001);
        
        if isempty(find1)==0
            Con_MAT_per(find1,i)=1;
            Con_MAT_per(i,find1)=1;
        end
        if isempty(find2)==0
            Con_MAT_per(find2,i)=1;
            Con_MAT_per(i,find2)=1;
        end
        
       
    end
end

for i=1:L
    I3=find((vertices(:,1)-Vcell(i,1)-(-p(1)+q(1))).^2+(vertices(:,2)-Vcell(i,2)-(-p(2)+q(2))).^2<0.0000001);
    L2=length(I3);
    for j=1:L2
        if mod(I3(j),4)==0
            v1=I3(j)-1; % the vs correspond to indices we join in the connectivity matrix
            
            v2=I3(j)-3;
            
        elseif mod(I3(j),4)==1
            v1=I3(j)+1;
            
            v2=I3(j)+3;
            
        elseif mod(I3(j),4)==2
            v1=I3(j)-1;
            v2=I3(j)+1;
            
        else
            
            v1=I3(j)-1;
            v2=I3(j)+1;
        end
        find1=find((vertices(v1,1)-Vcell(:,1)).^2+(vertices(v1,2)-Vcell(:,2)).^2<0.0000001);
        find2=find((vertices(v2,1)-Vcell(:,1)).^2+(vertices(v2,2)-Vcell(:,2)).^2<0.0000001);
        
        if isempty(find1)==0
            Con_MAT_per(find1,i)=1;
            Con_MAT_per(i,find1)=1;
        end
        if isempty(find2)==0
            Con_MAT_per(find2,i)=1;
            Con_MAT_per(i,find2)=1;
        end
        
       
    end
end



A = -Con_MAT_per;
A = A-spdiags(squeeze(sum(A,2)),0,size(A,1),size(A,2));
d = symrcm(A);
A = A(d,d);
end

