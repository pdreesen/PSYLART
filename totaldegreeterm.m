function tdeg = totaldegreeterm(term, vars),
% totaldegreeterm(term, vars)
% computes the total degree of the term 'term' wrt. 'vars'

if length(term) > 1,
    warning( 'several terms: considering only first term' );
    term = term(1);
end

if nargin < 2, 
    vars = symvar(term);
end



mdeg = multidegree(term, vars);

tdeg = sum(mdeg);

end

