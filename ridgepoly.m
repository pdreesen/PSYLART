function [coeffs,A,b,res,monsbasis]=ridgepoly(X,y,degree,gamma),
% Regularized polynomial regression
%
% [coeffs,A,b,res,monsbasis]=ridgepoly(X,y,degree,gamma)
%
% INPUTS
%   X: matrix with input data of size N x nvar 
%   y: vector with output data of size N x 1
%   degree: desired max degree of polynomial function
%   gamma: regularization constant / ridge parameter
%
% OUTPUTS 
%   coeffs: coefficients of fitted polynomial (ordered same as in
%   monsbasis)
%   A: matrix N x nbmons containing the evaluations of the input data X
%   b: vector N x 1 containing output data y
%   res: normwise residual A*coeffs-b
%   monsbasis: basis of monomials (indexing columns of A and coeffs
%
% EXAMPLE
%   X=[[1 2 4 4 12 2 4]' [2 3 4 2 4 1 5]'];
%   Y=[1 1 1 2 4 2 5]';
%   
%   deg=2;
%   gam=3;
%
%   [co,A,b,res,monbas]=ridgepoly(X,Y,deg,gam);
%   
%   [co monbas]
%
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   April 2012
%

nvar = size(X,2);
N = size(X,1);
nbmons = nb_mons_full(nvar,degree);
monsbasis = generate_mons_full(nvar,degree);

A=ones(N,nbmons);

for k = 2:nbmons, % eerste monomiaal is altijd (0,...,0)
    for i = 1:N,
        for j = 1:nvar,
            A(i,k) = A(i,k)*X(i,j)^monsbasis(k,j);
        end
               
    end
    
end

b=y;

Awo1=A(:,2:end);

coeffs=ridge(b,Awo1,gamma,0);

br=A*coeffs;

res = norm(br-b);

end

function nb = nb_mons_full(n,d),
% Return the number of monomials in a full monomial basis.
%  
% SIGNATURE
% nb=nb_mons_full(nvar,d);
%
% DESCRIPTION
% Returns the number of monomials in a full monomial basis of nvar unknowns 
% and total degree d (full: degrees 0:d)
%
% INPUTS
%    nvar   =   number of variables in input system
%    d      =   maximal total degree of full monomial basis
%
% OUTPUTS
%    nb     =   number of monomials in full monomial basis of nvar variables and degree d
%
% EXAMPLE
%
% CALLS
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010
%

nb=0;
for i=0:d,
	nb = nb + nb_mons_partial(n,i);
end

end

function fullBase = generate_mons_full(n,d)
% Generate a full basis of monomials.
%  
% SIGNATURE
% fullBase = generate_mons_full(nvar,d)
%
% DESCRIPTION
% Generates full basis (degrees 0:d) of monomials as a matrix;
% Each row of base refers to a n-tuple of exponents of a monomial 
% where each column corresponds to a variable.
% 
% INPUTS
%    nvar       =    number of variables
%    d          =    desired degree
%
% OUTPUTS
%    fullBase   =    full basis of monomials 
%
% EXAMPLE
%
% CALLS
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010
%

fullBase = zeros(1,n);

for i = 1 : d,
    fullBase = [fullBase ; generate_mons_partial(n,i)];
end

end

function base = generate_mons_partial(n,d)
% Generate a partial basis of monomials.
%  
% SIGNATURE
% base = generate_mons_partial(nvar,d)
%
% DESCRIPTION
% Generate a partial basis (degree 'd'; 'n' variables) of monomials as a 
% matrix. Each row of 'base' refers to a n-tuple of exponents of a monomial 
% where each column corresponds to a variable.
% 
% INPUTS
%    nvar       =    number of variables
%    d          =    desired degree
%
% OUTPUTS
%    base       =    partial basis of monomials 
%
% EXAMPLE
%
% CALLS
% generate_mons_partial
%
% AUTHOR
%    Philippe Dreesen (philippe.dreesen@gmail.com)
%    June 2010
%
% TODO: preallocate base and work with index of writeatrow or something

if n <= 1
    base = d;
else
    base = [];
    for i = d :-1: 0
        temp = generate_mons_partial(n-1,d-i);
        base = [base; i*ones(size(temp,1),1) generate_mons_partial(n-1,d-i)];
    end
end

end

function nb=nb_mons_partial(n,d),
% Return the number of monomials in a partial monomial basis.
%  
% SIGNATURE
% nb=nb_mons_partial(nvar,d);
%
% DESCRIPTION
% Returns the number of monomials in a partial monomial basis of nvar unknowns 
% and total degree d (partial: only degree d) 
%
% INPUTS
%    nvar   =   number of variables in input system
%    d      =   degree of partial monomial basis
%
% OUTPUTS
%    nb     =   number of monomials in partial monomial basis of nvar variables and degree d
%
% EXAMPLE
%
% CALLS
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010
%

nb = (factorial(d+n-1)/factorial(d))/factorial(n-1);

end