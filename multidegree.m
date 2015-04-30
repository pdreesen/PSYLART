function mdeg = multidegree(term, vars),
% mdeg = multidegree(term, [vars]);
% computes multidegree of (symbolic) term 'term' wrt. variables 'vars'


if length(term) > 1,
    warning( 'several terms: considering only first term' );
    term = term(1);
end

if nargin < 2, 
    vars = symvar(term);
end

mdeg = zeros(length(vars),1);

for i = 1:length(vars),
    mdeg(i) = degreei(term, vars(i));
end

end

