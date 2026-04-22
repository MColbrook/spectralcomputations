 function C = set_sum_lowmem_B(A,B,tol)

% function C = set_sum_lowmem(A,B)
%
% A, B are sets of intervals, of dimension (# intervals)-by-2
% C is returned as the sum of the intervals in A and B.
% C should not contain any overlapping intervals.
% The intervals of C are sorted from left to right.
% 
% This "low memory" version avoids creating storage for 
% length(A)*length(B) entries; it processes entries of A 
% one at a time, to minimize the creation of excessive 
% storage.

 lA = size(A,1);

 C=[];
 for j=1:lA
    CC = [A(j,1)+B(:,1) A(j,2)+B(:,2); C];        % add a new batch of intervals
    C=consolidate_int(CC,tol);
    
 end
