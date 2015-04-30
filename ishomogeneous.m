function res  = ishomogeneous(IN),
% Return whether a given polynomial system is homogeneous.
%
% SIGNATURE
% res = ishomogeneous(polyorig);
%
% DESCRIPTION
% Returns 1 if a system (or single equation) of polynomials 'polyorig' is 
% homogeneous. Returns 0 if the system is not homogeneous. 
%
% INPUTS
%    polyorig =   input system of polynomials (can be one equation given in
%                 matrix format)
%
% OUTPUTS
%    res      =   1 if the system is homogeneous, 0 otherwise
%
% EXAMPLE
%    ishomogeneous([1 2 0 1; 1 1 1 1]) returns 1
%    ishomogeneous([1 1 0 0; 1 2 0 1]) returns 0
%
% CALLS
%    get_info 
%
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   March 2011
%

% reformat input if necessary
if iscell(IN), 
    polyorig=IN;
elseif (~isscalar(IN) & ismatrix(IN)),
    polyorig{1} = IN;
else 
    error('unknown input format'); 
end

% get info
[neq, nvar, degrees, dmin] = get_info(polyorig);

% run over all equations and see whether for each term (=row) the degree
% is equal to degrees(eqni)
res=1;
for eqni = 1:neq,	
	for termi = 1:size(polyorig{eqni},1),
		degreedifference = degrees(eqni) - sum(polyorig{eqni}(termi,2:end));
		if degreedifference ~= 0,
            res =  0;
        end
	end
end


end % function
