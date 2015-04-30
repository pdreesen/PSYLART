function [rowsidc]=find_linindrows(Z,tol)
% Return the (row) indices of the linearly independent rows in a given matrix.
%  
% SIGNATURE
% [rowsidc]=find_linindrows(Z,tol)
%
% DESCRIPTION
% Returns the linearly independent rows of a given matrix, starting from
% the top.
% 
% INPUTS
%    Z          =    given matrix
%    tol        =    tolerance to be used in computing ranks
%
% OUTPUTS
%    rowsidc    =    the rows (starting from the top of Z) which are
%                    linearly independent; the length of rowsidc
%                    corresponds to the rank of Z
%
% EXAMPLE
%  
% CALLS
% 
% AUTHOR
%     Philippe Dreesen (philippe.dreesen@gmail.com)
%     September 2011
%

if nargin < 2, tol = 1e-10; end

rankslist=zeros(size(Z,1),1);
ranksdiff=zeros(size(Z,1),1);

for i = 1:size(Z,1),
	rankslist(i) = rank(Z(1:i,:),tol);
    
	if (i==1), 
        ranksdiff(1) = rankslist(1);
    elseif (i>1), 
        ranksdiff(i) = rankslist(i) - rankslist(i-1);
	end
end	
rowsidc=find(ranksdiff);

end