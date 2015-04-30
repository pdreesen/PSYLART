function polysys=polyorig_to_polysys(polyorig),
% Convert polynomial system in polyorig format to polysys format.
%  
% SIGNATURE
% polysys=polyorig_to_polysys(polyorig)
%
% DESCRIPTION
% Converts polynomial system in polyorig format to polysys format 
% 
% INPUTS
%    polyorig   =    system of polynomials in polyorig format (a neq x 1 cell 
%                    array with corresponding to each input equation a matrix 
%                    with for each term a row of the form [coeff exps], where 
%                    coeff is the coefficient and exps is a representation of 
%                    the monomials (e.g., 1 1 for xy)
%
% OUTPUTS
%    polysys    =    system of polynomials in polysys format (a neq x 2 cell
%                    array with in the {i,1} entry a vector containing
%                    the coefficients of input equation i, and in the {i,2} 
%                    entry a no_terms x nvar matrix containing the exponents 
%                    of the monomials for all terms in the input equation i. 
%
% EXAMPLE
%
% CALLS
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010
%


neqs = size(polyorig,2);
polysys = cell(neqs,2);

for i=1:neqs,
	polysys{i,1} = polyorig{i}(:,1)';
	polysys{i,2} = polyorig{i}(:,2:end);
end


end
