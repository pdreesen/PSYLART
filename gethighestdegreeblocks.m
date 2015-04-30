function RES=gethighestdegreeblocks(Md,n,d,degs),
% Return the highest degree blocks of a given Macaulay matrix.
%  
% SIGNATURE
% RES=gethighestdegreeblocks(Md,nvar,d,degs)
%
% DESCRIPTION
% Returns the highest degree blocks of a given Md matrix M
%
% INPUTS
%    Md         =    coefficient matrix
%    nvar       =    number of variables
%    d          =    degree of coefficient matrix Md
%    degs       =    list containing the desired degree blocks 
%                    (e.g., if d=8, degs=5:8 will return the blocks
%                    corresponding to degrees 5 up to 8)
%
% OUTPUTS
%    RES        =    degs-degree blocks from M 
%
% EXAMPLE
%
% CALLS
% generate_mons_partial
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010
%


% function RES= gethighestdegreeblocks(Md,n,d,degs);
% returns the highest degree blocks of a given Md matrix M
% Md: computed Md matrix
% n: number of variables in original system from which Md is generated
% d: degree of Md matrix
% degs: list containing the desired degree blocks (e.g., if d=8, degs=5:8 will return the blocks
% corresponding to degrees 5 up to 8)

numberofcolsfromend=0;
for i=1:length(degs),
	numberofcolsfromend=numberofcolsfromend + size(generate_mons_partial(n,degs(i)),1);
	size(generate_mons_partial(n,i),1);
end
RES = Md(:,end-numberofcolsfromend+1:end);
end
