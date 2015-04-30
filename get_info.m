function [neq, nvar, degrees, dmin, coeffs, expons, bezout] = get_info(polyorig)
% Extract information from a given polynomial system.
%  
% SIGNATURE
% [neq, nvar, degrees, dmin, coeffs, expons, bezout] = get_info(polyorig)
%
% DESCRIPTION
% Extract information from the polynomial system polyorig
% 
% INPUTS
%    polyorig    =   input system of polynomials
%
% OUTPUTS
%   neq          =   number of equations
%   nvar         =   number of variables
%   degrees      =   degrees of input equations
%   dmin         =   minimal degree of which M must be (that is, dmin=max(sum(degrees)))
%   coeffs       =   cell containing the coefficients of the system for each equation
%   expons       =   cell containing the exponents of the system for each equation
%   bezout       =   bezout number of system of polynomial equations
%
% EXAMPLE
%
% CALLS
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010
%

% ensure input is cell
if ~iscell(polyorig),
    pmat=polyorig;
    clear polyorig;
    polyorig{1}=pmat;
end
    
% number of equations
neq=length(polyorig);

% number of variables
nvar=size(polyorig{1},2)-1;

% degrees, coefficients, exponents
degrees = zeros(neq,1);
coeffs = cell(neq,1);
expons = cell(neq,1);
for i = 1:neq,
	    degrees(i) = max(sum(polyorig{i}(:,2:end),2));
		coeffs{i} = [polyorig{i}(:,1)];
		expons{i} = [polyorig{i}(:,2:end)]; 
end

% minimal degree (=max total degree of original system)
dmin=max(degrees);

% bezout count (=multiplication of degrees) 
bezout = 1;
for i=1:neq, 
    bezout=bezout * degrees(i);
end

end
