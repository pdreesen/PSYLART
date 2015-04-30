function out=partderiv(in,degrees),
% Derive a collection of terms with respect to several variables and given degrees.
% SIGNATURE
% out=partderiv(in,degrees);
%
% DESCRIPTION
% Derives a collection of terms (i.e., a matrix where each row
% corresponds to one term) with respect to several variables and 
% given orders.
% 
% INPUTS 
%
% in        =   matrix containing rows of the form [coeff exps], 
%               where 'coeff' is a scalar depicting the coefficient 
%               of the monomial, and 'exps' are the powers of the
%               variables
%
% degrees   =   the derivation orders for the variables 
%
% OUTPUT
% out       =   the rows after computing the derivatives
%
% EXAMPLE
% If one calls 
%               out = partderiv([2 2 1;1 3 2;1 1 1],[1 2])
% the terms in the first argument represent 2 x^2 y, x^3 y^2 and x y,
% and deriving these with respect to the orders [1 2], (i.e., first
% order derivative in x, second order derivative in y), this yields
%
%   out =
%   
%        0     0     0
%        6     2     0
%        0     0     0
%
% CALLS
%       deriveallrows1varhighdegree
% AUTHOR
%       Philippe Dreesen (philippe.dreesen@gmail.com)
%       KULeuven, ESAT/SCD
%       July 2010
% TODO:
% normalization of partderiv: divide by (a1!*a2!*...*an!) where
% (a1,a2,...,an) is the derivation orders


out=in;
nvars=size(in,2)-1;
for i =1:nvars,
	out = deriveallrows1varhighdegree(out,i,degrees(i));	
end

out(:,1) = out(:,1)/prod(factorial(degrees));

end
















function out=deriveallrows1varhighdegree(in,var,deg),
% SIGNATURE
% out=deriveallrows1varhighdegree(in,var,deg)
%
% DESCRIPTION
% Derives a collection of terms (i.e., a matrix where each row
% corresponds to one term) with respect to one variable (i.e., 
% also only one order) and arbitrary (>= 1) order. Also see d
% derive1row1var for further info.
% 
% INPUTS 
% in        =   matrix containing rows of the form [coeff exps], 
%               where 'coeff' is a scalar depicting the coefficient 
%               of the monomial, and 'exps' are the powers of the
%               variables
%
% var       =   the the i-th variable with respect to which the
%               derivative is to be computed
%
% deg       =   the degree to which the terms have to be derived
%               (note that this is an integer, and corresponds to 
%               derivation with respect to the variable 'var')
%
% OUTPUT
% out       =   the terms (rows) after computing the derivative
%
% EXAMPLE
% If one calls 
%               out = deriveallrows1varhighdegree([4 2 1;5 1 2],2,2)
% the terms in the first argument represent 4 x^2 y and 5 x y^2, 
% and deriving this with respect to the second variable (var=2), 
% being y, yields        
%    out =
%    
%         0     0     0
%        10     1     0
%
%
% CALLS
%       deriveallrows1var
%
% AUTHOR
%       Philippe Dreesen (philippe.dreesen@gmail.com)
%       KULeuven, ESAT/SCD
%       July 2010
%

out=in;

for i=1:deg,
	out = deriveallrows1var(out,var); 
end

end




















function out=deriveallrows1var(in,var),

% SIGNATURE
% out=deriveallrows1var(in,var);
%
% DESCRIPTION
% Derives a collection of terms (i.e., a matrix where each row
% corresponds to one term) with respect to one variable (i.e., 
% also only one order). Also see derive1row1var for further info.
% 
% INPUTS 
% in        =   matrix containing rows of the form [coeff exps], 
%               where 'coeff' is a scalar depicting the coefficient 
%               of the monomial, and 'exps' are the powers of the
%               variables
%
% var       =   the the i-th variable with respect to which the
%               derivative is to be computed
%
% OUTPUT
% out       =   the rows after computing the derivative
%
% EXAMPLE
% If one calls 
%               out = derive1row1var([4 2 1;5 1 2],1)
% the terms in the first argument represent 4 x^2 y and 5 x y^2, 
% and deriving this with respect to the first variable (var=1), 
% being x, yields        
%    out =
%    
%         8     1     1
%         5     0     2
%
%
% CALLS
%       derive1row1var
%
% AUTHOR
%       Philippe Dreesen (philippe.dreesen@gmail.com)
%       KULeuven, ESAT/SCD
%       July 2010
%

out=zeros(size(in));

for i=1:size(in,1),
	out(i,:)=derive1row1var(in(i,:),var);
end

end


















function out=derive1row1var(in,var),

% SIGNATURE
% out=derive1row1var(in,var);
%
% DESCRIPTION
% Derives a term (i.e., one row only) with respect to one variable (i.e., also only one order).
% 
% INPUTS 
%
% in        =   one term (row vector) of the form [coeff exps], where 
%               'coeff' is a scalar depicting the coefficient of the 
%               monomial, and 'exps' are the powers of the variables
%
% var       =   the the i-th variable with respect to which the
%               derivative is to be computed
%
% OUTPUT
% out       =   the term after computing the derivative
%
% EXAMPLE
% If one calls 
%               out = derive1row1var([4 2 1],1)
% the term in the first argument represents 4 x^2 y, and deriving 
% this with respect to the first variable (var=1), being x, this 
% yields        
%               out=[8 1 1]
%
% CALLS
%
% AUTHOR
%       Philippe Dreesen (philippe.dreesen@gmail.com)
%       KULeuven, ESAT/SCD
%       July 2010
%

out=zeros(1,size(in,2));

if in(var+1)~=0, 
	out = in;
	out(1)=in(1)*in(var+1);
	out(var+1)=in(var+1)-1;
end

end
