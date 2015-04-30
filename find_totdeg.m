function tdeg=find_totdeg(polyvec,nvar),
% Return the total degree of a polynomial.
%  
% SIGNATURE
% tdeg=find_totdeg(polyvec,nvar)
%
% DESCRIPTION
% Returns the total degree of the polynomial polyvec
%
% INPUTS
%    polyvec    =    input polynomial equation
%    nvar       =    number of variables
%
% OUTPUTS
%    tdeg       =    total degree of equation polyvec
%
% EXAMPLE
%
% >> find_totdeg([1,2,3,4,5,6],2)
% 
% ans =
% 
%      2
%
% CALLS
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010

 
lpolyv=length(polyvec);
tdeg = 0;
while nb_mons_full(nvar,tdeg)<lpolyv,
	tdeg=tdeg+1;
end


end
