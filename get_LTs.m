function polyOUT=get_LTs(polyorig),
% Return the highest degree terms of a given system of polynomial
% equations.
%
% SIGNATURE
% polyOUT = get_LTs(polyorig)
%
%
% DESCRIPTION
% Return a system consisting of only the highest degree terms of a given
% system of polynomial equations. 
%
%
% INPUTS
%    polyIN   =  input system of max total degree deg
%
%
% OUTPUTS
%    polyOUT  =  output system consisting of only the highest degree terms
%
%
% EXAMPLE
% >> polyorig{1} = [1 2 0 0;1 0 0 5;3 1 1 1;5 1 0 4];
% >> polyOUT=get_LTs(polyorig);
% 
%
% CALLS
% get_info, nb_mons_full, vec_to_polyorig, polyorig_to_vec
%
%
% AUTHOR
%       Philippe Dreesen (philippe.dreesen@gmail.com)
%       KU Leuven, ESAT/SCD
%       May 2012
%


[neq, nvar, degrees, dmin, coeffs, expons, bezout]=get_info(polyorig);

polyOUT=cell(neq,1);

for i = 1:neq,
   pvec=polyorig_to_vec(polyorig{i});
   pvec(1:nb_mons_full(nvar,degrees(i)-1)) = 0;
   polyOUT{i}=vec_to_polyorig(pvec,nvar);
end

end