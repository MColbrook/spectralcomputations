function [C] = consolidate_int(CC,tol)
    [~,indx] = sort(CC(:,1),'ascend');                     % sort by left end
    CC = CC(indx,:);
    C=CC;


    k = 1;
    lCC = length(CC(:,1));
    cpos = 0;
    while (k <= lCC)                              % consolidate intervals
       % left  = CC(k,1)-tol;
       right = CC(k,2);
       jj = k+1;
       while ( (jj <= lCC) && (CC(jj,1) <= right+tol) )  % find all overlapping intervals
          right = max(right,CC(jj,2)); 
          jj = jj+1;
       end
       cpos = cpos+1;
       C(cpos,1) = CC(k,1);                          % this indexing is slightly faster
       C(cpos,2) = right;
       k = jj;
    end
    C = C(1:cpos,:);  
end