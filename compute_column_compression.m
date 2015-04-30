function Acc=compute_column_compression(A,comprank,rowind)
% Compute a column compressed version of a given matrix.
%
% SIGNATURE
% Acc=compute_column_compression(A,naff,rowind);
% 
% DESCRIPTION
% Compute the column compression of a given matrix A given the number of 
% columns (estimation of rank) to retain and an index to where the 
% A(1:rowind,:) part is assumed to span the estimated rank.
%  
% INPUTS
%    A         =   matrix
%    comprank  =   number of columns to retain which are assumed to carry
%                  rank.
%    rowind    =   row index such that A(1:rowind,:) is assumed to have 
%                  numerical rank comprank
% 
% OUTPUTS
%    Acc       =   column compressed version of A
%
% EXAMPLE
%    TODO
%
% CALLS
% svd
%
% AUTHOR
%    Philippe Dreesen (philippe.dreesen@gmail.com)
%    KULeuven, ESAT/SCD
%    June 2012
%

%column compression on A
[~,~,W]=svd(A(1:rowind,:));
Acc = A*W;
Acc=Acc(:,1:comprank);


end

