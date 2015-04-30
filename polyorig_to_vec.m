function pvec=polyorig_to_vec(polyorig,deg),
% Convert polyorig equation into vector format.
%
% SIGNATURE
% pvec=polyorig_to_vec(polyorig,deg)
% 
% DESCRIPTION
% Convert a given polynomial equation in polyorig format into vector
% format. 
%
% INPUTS
%    polyorig =   input polynomial equation
%    deg      =   max degree of vector representation (default is max
%                 degree of the input equation)
% 
% OUTPUTS
%    pve      =   vector representation of input polynomial
%
% EXAMPLE
% >> pvec=polyorig_to_vec([1 0 2;2 1 1;3 2 0;4 0 1; 5 1 0; 6 0 0])
% 
% pvec = 
%      6
%      5
%      4
%      3
%      2
%      1
%
% CALLS
%       get_info, nb_mons_full, fetr
%
% AUTHOR
%       Philippe Dreesen (philippe.dreesen@gmail.com)
%       KU Leuven, ESAT/SCD
%       May 2012
%

% ensure input is a matrix
% if input is cell, take contents, which should be a matrix
if iscell(polyorig),
    pmat=polyorig{1};
elseif isnumeric(polyorig)
    pmat=polyorig;
else
    error('input format polynomial is wrong!');
end

% obtain info of equation
[neq, nvar, degrees, dmin, coeffs, expons, bezout] = get_info(pmat);

if nargin<2,
    pvecsize=nb_mons_full(nvar,dmin);
else
    pvecsize=nb_mons_full(nvar,deg);
end

pvec=zeros(pvecsize,1);

for i = 1:size(pmat,1),
    E=pmat(i,2:end);
    idx=fetr(E);
    pvec(idx)=pvec(idx)+pmat(i,1);    
end


end