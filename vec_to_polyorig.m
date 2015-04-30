function polyorig = vec_to_polyorig(polyvec,nvar,tol),
% Convert a polynomial from vector format to polyorig format.
%  
% SIGNATURE
% polyorig = vec_to_polyorig(polyvec,nvar,tol)
%
% DESCRIPTION
% Converts a polynomial from vector format to polyorig format, given the
% number of variables used
%
% INPUTS
%    polyvec    =    input polynomial equation
%    nvar       =    number of variables
%    tol        =    tolerance used for detecting zero elements
%
% OUTPUTS
%    polyorig   =    polyorig representation of the polyvec polynomial
%
% EXAMPLE
%
% vec_to_polyorig([1,2,3,4,5,6],2)
% 
% ans =
% 
%      1     0     0
%      2     1     0
%      3     0     1
%      4     2     0
%      5     1     1
%      6     0     2
%
%
% CALLS
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010/May 2011/July 2011

% default tolerance
if (nargin<3), tol=1e-8;end

polyvec=polyvec(:)';

% detecting zero elements
polyvec(abs(polyvec)<tol)=0;

% main part
tdeg = find_totdeg(polyvec,nvar);

monsfull = generate_mons_full(nvar,tdeg);

polyorig = zeros(1,nvar+1);
writeatrow = 1;

for i=1:length(polyvec),
	if polyvec(i)~=0,
		polyorig(writeatrow,1)=polyvec(i);
		polyorig(writeatrow,2:end)=monsfull(i,:);
		writeatrow=writeatrow+1;
	end
end

end
