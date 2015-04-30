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


% fullBase = generate_mons_full(n,d)
% Generates (full) basis of monomials of degree d in n variables.


fullBase = zeros(1,n);

for i = 1 : d,
    fullBase = [fullBase ; generate_mons_partial(n,i)];
end

end
