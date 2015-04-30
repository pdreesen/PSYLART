function K=build_canonical_kernel(nvar,d,sols),
% Build canonical kernel given the solutions.
%  
% SIGNATURE
% K=build_canonical_kernel(nvar,d,sols);
%
% DESCRIPTION
% Build a canonical kernel 'K' having columns consisting of 
% the canonical monomial form evaluated at the roots given 
% by 'sols' with total degree 'd'
% 
% INPUTS
%    nvar    =    number of variables
%    d       =    requested degree of monomial basis
%    sols    =    cell containing solution n-tuplets in vectors
%
% OUTPUTS
%    K       =    canonical kernel evaluated at the given roots
%
% EXAMPLE
% For the case nvar = 2, d=3, sols = {(1,2),(3,4)}, this yields 
%        [ 1  1 ]  1
%        [ 1  3 ]  x
%        [ 2  4 ]  y
%        [ 1  9 ]  x^2
%    K = [ 2 12 ]  x y
%        [ 4 16 ]  y^2
%        [ 1 27 ]  x^3
%        [ 2 36 ]  x^2 y
%        [ 4 48 ]  x y^2
%        [ 8 64 ]  y^3
%
% CALLS
% nb_mons_full, generate_mons_full
% 
% AUTHOR
%    Philippe Dreesen (philippe.dreesen@gmail.com)
%    KULeuven, ESAT/SCD
%    June 2010
%

rowsK = nb_mons_full(nvar,d);
colsK = length(sols);

K = ones(rowsK,colsK);

basiske = generate_mons_full(nvar,d);

for colcol = 1:colsK,
	for rowrow=1:rowsK,
		for vari = 1:nvar,
			K(rowrow,colcol) = K(rowrow,colcol) * sols{colcol}(vari).^basiske(rowrow,vari);
		end
	end
end

end % function



