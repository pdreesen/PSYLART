function pv = polysym_to_polyvec(polysym,vars),
% convert symbolic (MATLAB MuPAD) polynomial to polynomial vector (degree
% negative lexicographic order)

if nargin<2,
    vars = symvar(polysym);
end

[c,t] = coeffs(polysym,vars);
c = double(c);

tdeg = totaldegree(polysym);

mons = generate_mons_full( length( vars ), tdeg );
pv = zeros( nb_mons_full( length( vars ), tdeg ), 1 );

for ti = 1: length(c),
    monti = multidegree( t(ti), vars );
    pv( find_mon_index( monti ,mons ) ) = c(ti);
end

end