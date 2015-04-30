function [cor,s] = corank(MATRIX,tol,lr)
% Compute the corank of a matrix.
%
% SIGNATURE
% cor = corank(M, tol, lr)
% 
% DESCRIPTION
% Computes the (left or right) corank of matrix MATRIX, using threshold 
% tol in computing the rank (using an SVD).
%  
% INPUTS
%    M        =   given matrix
%    tol      =   tolerance used in determining numerical rank (in SVD)
%    lr       =   flag 'l' or 'r' (default) denoting whether the left ('l')
%                 or right ('r') corank has to be computed
% 
% OUTPUTS
%    cor      =   corank of given matrix M
%    s        =   vector containing singular values of M 
%
% EXAMPLE
% >> A = [ 1 0 -2 0 4; 0 1 6 0 1; 1 0 1 0 -4]
% 
% A =
% 
%      1     0    -2     0     4
%      0     1     6     0     1
%      1     0     1     0    -4
% 
% >> corank(A)
% 
% ans =
% 
%      2
%
% CALLS
%
% SEE ALSO
% svd  
% 
% AUTHOR
%    Philippe Dreesen (philippe.dreesen@gmail.com)
%    KULeuven, ESAT/SCD
%    June 2010
%

s = svd(MATRIX);

% by default take right corank (=dim of right nullspace)
if nargin <3, 
    lr = 'r';
end

% check if tolerance is set, if not, set!
if (nargin < 2) || isempty(tol), 
	tol = max(size(MATRIX)) * eps(max(s));
end

% compute rank
r = sum(s > tol);

% compute `column dimension' (can be number of cols in case of right corank, 
% or the number of rows in case of left corank)
if lr == 'r',
    coldim = size(MATRIX,2);
elseif lr == 'l',
    coldim = size(MATRIX,1);
else
    error('unknown stuffs');
end

% compute corank as column dimension minus rank
cor = coldim -r ;
end
