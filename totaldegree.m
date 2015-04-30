function tdeg = totaldegree(F,vars),
% tdeg = totaldegree(F, [vars]);
% computes total degree of a polynomial (over all terms) with respect 
% to the variables 'vars'

F=sym(F);

if nargin < 2,
    vars = symvar(F);
end

tdeg = 0;
[c,t]=coeffs( F );
% c = double(c);

for ti = 1:length(c),
    tdeg = max(tdeg,totaldegreeterm(t(ti),vars));
end



end