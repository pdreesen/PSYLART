function out = evalfunctional(func,root),
% Evaluate a collection of terms in a given root.
%
% SIGNATURE
% out=evalfunctional(func,root);
%
% DESCRIPTION
% Evaluates a collection of terms (i.e., a matrix where each row
% corresponds to one term of the form [coeff exps]) in a given
% root (of the form [x1* x2* ... xn*].
% 
% INPUTS 
%
% func      =   matrix containing rows of the form [coeff exps], 
%               where 'coeff' is a scalar depicting the coefficient 
%               of the monomial, and 'exps' are the powers of the
%               variables
%
% root      =   the given root in which the functional is to be 
%               evaluated 
%
% OUTPUT
% out       =   the rows after evaluating the given root 
%
% EXAMPLE
% If one calls 
%               out = evalfunctional([2 2 1;1 3 2;1 1 1],[1 2])
% the terms in the first argument represent 2 x^2 y, x^3 y^2 and x y,
% and evaluating these terms at the given root [1 2], this yields
%
%   out =
%
%        4
%        4
%        2
%
% CALLS
% evalfunctionalrow       
%
% AUTHOR
%    Philippe Dreesen (philippe.dreesen@gmail.com)
%    KULeuven, ESAT/SCD
%    July 2010
%

for i = 1:size(func,1),
	out(i,:) = evalfunctionalrow(func(i,:),root);
end


end







function outrow = evalfunctionalrow(term,root),
% evaluates one row only

outrow = term(1);
for j=1:length(root),
	outrow = outrow*root(j).^term(1+j);
end

end



