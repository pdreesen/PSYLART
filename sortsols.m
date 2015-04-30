function solssorted=sortsols(sols),
% Sort the roots according to their modulus
%  
% SIGNATURE
% solssorted = sortsols(sols);
%
% DESCRIPTION
% Sorts the solutions in sols according to increasing modulus.
% 
% INPUTS
%    sols         =    solutions matrix of size nbsols x nbunknwns
%
% OUTPUTS
%    solssorted   =    sorted matrix of solutions
%
% EXAMPLE
% 
% CALLS
% 
% AUTHOR
%     Philippe Dreesen (philippe.dreesen@gmail.com)
%     March 2012
%

nbsols = size(sols,1);
nbunks = size(sols,2);

normsol=zeros(nbsols,1);

for i = 1:nbsols,
   normsol(i) = norm(sols(i,:),2);
end

[~,order]=sort(normsol);

solssorted=sols(order,:);

end

