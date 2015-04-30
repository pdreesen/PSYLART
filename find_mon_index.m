function idx=find_mon_index(newmon,mons)
% Find the index of a monomial in a list of monomials.
%
% SIGNATURE
% idx=find_mon_index(newmon,mons)
% 
% DESCRIPTION
% Find the index of a monomial 'newmon' in a list of monomials 'mons'
%
% INPUTS
%    newmon   =   monomial contained in mons at unknown row
%    mons     =   given list of monomials (each row corresponds
%                 to a monomial) 
% 
% OUTPUTS
%    idx      =   row number at which newmon occurrs in mons
%
% EXAMPLE
%       TODO
%
% CALLS
%       TODO
%
% AUTHOR
%       Philippe Dreesen (philippe.dreesen@gmail.com)
%       KULeuven, ESAT/SCD
%       June 2010
%
% TODO:
% use fetr, frte

% idx=find_mon_index(newmon, mons)
% Find the position (i.e., index or row number) where the monomial newmon occurs in the basis of monomials mons

idx = 0;
for posi=1:size(mons,1),
	if (isequal(newmon,mons(posi,:))),
		idx=posi;
	end
end

if idx == 0, 
    disp('something went wrong; index of monomial could not be found');
end

end
