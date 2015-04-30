function polyorig=polysys_to_polyorig(polysys),
% Convert a given polynomial system in polysys format to polyorig format.
%  
% SIGNATURE
% polyorig=polysys_to_polyorig(polysys)
%
% DESCRIPTION
% Converts polynomial system in polysys format to polyorig format 
% 
% INPUTS
%    polysys    =    system of polynomials in polysys format (a neq x 2 cell
%                    array with in the {i,1} entry a vector containing
%                    the coefficients of input equation i, and in the {i,2} 
%                    entry a nterms x nvar matrix containing the exponents 
%                    of the monomials for all terms in the input equation i. 
%
% OUTPUTS
%    polyorig   =    system of polynomials in polyorig format (a neq x 1 cell 
%                    array with corresponding to each input equation a matrix 
%                    with for each term a row of the form [coeff exps], where 
%                    coeff is the coefficient and exps is a representation of 
%                    the monomials (e.g., 1 1 for xy)
%
% EXAMPLE
%
% CALLS
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010
%


neqs = size(polysys,1);
polyorig = cell(1,neqs);

for i=1:neqs,
	polyorig{i} = [polysys{i,1}' polysys{i,2}]; 
end


end
