function polyorighom=homogenize_polyorig(polyorig),
% Homogenize a given polynomial system.
%  
% SIGNATURE
% polyorighom=homogenize_polyorig(polyorig);
%
% DESCRIPTION
% Homogenizes a system of polynomial equations: returns homogeneous 
% system which is built by multiplying each term with a newly introduced
% dummy variable raised to a power, such that every term in 
% each equation in the resulting system is of equal degree
% 
% INPUTS
%    polyorig   =    input system of polynomials
%
% OUTPUTS
%    polyorigh  =    system of homogeneous polynomials 
%
% EXAMPLE
% >> polyorig{1} = [1 1 2; -3 2 0; 1 1 0; -2 0 1]; 
% >> polyorig{2} = [2 1 2; -1 0 1; 1 1 0; -5 0 0];
% >> polyhom = homogenize_polyorig(polyorig);
% >> polyhom{:} 
% 
% ans =
% 
%      1     1     2     0
%     -3     2     0     1
%      1     1     0     2
%     -2     0     1     2
% 
% 
% ans =
% 
%      2     1     2     0
%     -1     0     1     2
%      1     1     0     2
%     -5     0     0     3
% 
% 
% CALLS
%    get_info 
%
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010
%

% get info from system
[neq, nvar, degrees, dmin, coeffs, expons, bezout]=get_info(polyorig);

% initialize homogeneous system
polyorighom = cell(neq,1);

for eqni = 1:length(polyorig),	
	% for each equation, loop over all terms and homogenize
	di = degrees(eqni);
	writeatrow = 1;
	for termi = 1:size(polyorig{eqni},1),
		degreedifference = di - sum(polyorig{eqni}(writeatrow,2:end));
		polyorighom{eqni}(writeatrow,:) = [polyorig{eqni}(writeatrow,:) degreedifference];
		writeatrow=writeatrow+1;
	end
end

end %function
