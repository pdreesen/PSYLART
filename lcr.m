function [cor,s] = lcr(A,tol),
% Compute the left corank of a matrix.
%
% SIGNATURE
% cor = lcr(A, tol)
% 
% DESCRIPTION
% Computes the left corank of matrix A, using threshold 
% tol in computing the rank (using an SVD).
%  
% INPUTS
%    A        =   given matrix
%    tol      =   tolerance used in determining numerical rank (in SVD)
% 
% OUTPUTS
%    cor      =   left corank of given matrix A
%    s        =   vector containing singular values of A
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
% >> lcr(A)
% 
% ans =
% 
%      0
%
% CALLS
% corank
%
% SEE ALSO
% svd  
% 
% AUTHOR
%    Philippe Dreesen (philippe.dreesen@gmail.com)
%    KULeuven, ESAT/SCD
%    June 2010
%

if nargin < 2, tol=[]; end

[cor,s] = corank(A,tol,'l')

end