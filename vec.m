function v = vec( x )

% VEC   Vectorize.
%    VEC(X), where X is a vector, matrix, or N-D array, returns a column vector
%    containing all of the elements of X; i.e., VEC(X)=X(:). 
%
% Copyright 2005-2013 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

v = reshape( x, numel( x ), 1 );


end


